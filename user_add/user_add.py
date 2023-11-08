import synapseclient as sc
import bridgeclient as bc
import pandas as pd
import random
import boto3
import json
import os
import re
from botocore.exceptions import ClientError

TESTING = False

def read_args():
    # for testing
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputTable")
    parser.add_argument("--outputTable")
    parser.add_argument("--bridgeUsername")
    parser.add_argument("--bridgePassword")
    parser.add_argument("--synapseUsername")
    parser.add_argument("--synapsePassword")
    parser.add_argument("--substudy")
    parser.add_argument("--supportEmail")
    args = parser.parse_args()
    return args


def get_env_var_credentials():
    credentials = {}
    credentials['inputTable'] = os.getenv('inputTable')
    credentials['outputTable'] = os.getenv('outputTable')
    credentials['synapseUsername'] = os.getenv('synapseUsername')
    credentials['synapsePassword'] = os.getenv('synapsePassword')
    credentials['bridgeUsername'] = os.getenv('bridgeUsername')
    credentials['bridgePassword'] = os.getenv('bridgePassword')
    credentials['substudy'] = os.getenv('substudy')
    credentials['supportEmail'] = os.getenv('supportEmail')
    return credentials


def delete_na_rows(syn, input_table):
    rows_to_delete = syn.tableQuery(
            "select * from {} where phone_number is null or guid is null".format(input_table))
    syn.delete(rows_to_delete)

def get_new_users(syn, input_table, output_table):
    input_table_df = syn.tableQuery(
            "select * from {}".format(input_table)).asDataFrame()
    for i, user in input_table_df.iterrows():
        if pd.isnull(user.phone_number) and pd.isnull(user.guid):
            delete_na_rows(syn, input_table)
            return ("Error: phone_number and guid was left blank", -1, -1, user.visit_date)
        elif pd.isnull(user.phone_number):
            delete_na_rows(syn, input_table)
            return ("Error: phone_number was left blank", -1, user.guid, user.visit_date)
        elif pd.isnull(user.guid):
            delete_na_rows(syn, input_table)
            return ("Error: guid was left blank", user.phone_number, -1, user.visit_date)
    phone_number_ints = input_table_df.phone_number.apply(get_phone_number_digits)
    input_table_df['phone_number'] = phone_number_ints
    input_table_df['phone_number'] = input_table_df.phone_number.astype(str)
    input_table_df['guid'] = input_table_df.guid.apply(str.strip)
    input_table_df = input_table_df.set_index(["phone_number", "guid"], drop=False)
    output_table_df = syn.tableQuery(
            "select phone_number, guid from {}".format(output_table)).asDataFrame()
    output_table_df['phone_number'] = output_table_df.phone_number.astype(str)
    output_table_df = output_table_df.set_index(["phone_number", "guid"], drop = False)
    new_numbers = set(input_table_df.index.values).difference(
            output_table_df.index.values)
    if len(new_numbers):
        return input_table_df.loc[list(new_numbers)]
    else:
        return input_table_df.drop(input_table_df.index)


def get_bridge_client(bridge_username, bridge_password, study="sage-mpower-2"):
    bridge = bc.bridgeConnector(
            email = bridge_username, password = bridge_password, study = study)
    return bridge


def get_participant_info(bridge, phone_number):
    participant_info = bridge.restPOST(
            "/v3/participants/search",
            {"phoneFilter": phone_number})
    return participant_info


def process_request(bridge, participant_info, phone_number, external_id,
                    substudy, support_email):
    if participant_info['total'] == 0:
        # create account
        try:
            # assign to random engagement group
            engagement_groups = []
            engagement_groups.append(random.choice(["gr_SC_DB", "gr_SC_CS"]))
            engagement_groups.append(random.choice(["gr_BR_AD", "gr_BR_II"]))
            engagement_groups.append(random.choice(["gr_ST_T", "gr_ST_F"]))
            engagement_groups.append(random.choice(["gr_DT_F", "gr_DT_T"]))
            bridge.restPOST(
                    "/v3/participants",
                    {"externalIds": {substudy: external_id},
                     "checkForConsent": False,
                     "phone": {"number": phone_number,
                               "regionCode": "US"},
                     "dataGroups": engagement_groups + ["clinical_consent", "show_3_cognitive"],
                     "sharingScope": "all_qualified_researchers"}) # assume US?
            return "Success: User account created"
        except Exception as e:
            return ("Error: Could not create user account. "
                    "Does your phone number have a US area code and/or "
                    "has the GUID already been assigned? "
                    "Console output: {0}".format(e))
    elif substudy not in participant_info['items'][0]['externalIds']:
        # add external_id and then assign to existing account
        user_id = participant_info['items'][0]['id']
        user_info = bridge.restGET("/v3/participants/{}".format(user_id))
        if "clinical_consent" not in user_info["dataGroups"]:
            user_info["dataGroups"] = user_info["dataGroups"] + ["clinical_consent"]
        if "show_3_cognitive" not in user_info["dataGroups"]:
            user_info["dataGroups"] = user_info["dataGroups"] + ["show_3_cognitive"]
        if not ("gr_SC_DB" in user_info["dataGroups"] or
                "gr_SC_CS" in user_info["dataGroups"]):
            user_info["dataGroups"] = user_info["dataGroups"] + \
                    [random.choice(["gr_SC_DB", "gr_SC_CS"])]
        if not ("gr_BR_AD" in user_info["dataGroups"] or
                "gr_BR_II" in user_info["dataGroups"]):
            user_info["dataGroups"] = user_info["dataGroups"] + \
                    [random.choice(["gr_BR_AD", "gr_BR_II"])]
        if not ("gr_ST_T" in user_info["dataGroups"] or
                "gr_ST_F" in user_info["dataGroups"]):
            user_info["dataGroups"] = user_info["dataGroups"] + \
                    [random.choice(["gr_ST_T", "gr_ST_F"])]
        if not ("gr_DT_F" in user_info["dataGroups"] or
                "gr_DT_T" in user_info["dataGroups"]):
            user_info["dataGroups"] = user_info["dataGroups"] + \
                    [random.choice(["gr_DT_F", "gr_DT_T"])]
        user_info["sharingScope"] = "all_qualified_researchers"
        try:
            bridge.restPOST("/v3/participants/{}".format(user_id), user_info)
        except Exception as e:
            return ("Error: Preexising user account found. "
                    "but data groups could not be updated. "
                    "Console output: {0}".format(e))
        if not "at-home-pd" in participant_info['items'][0]['externalIds']:
            try:
                bridge.restPOST("/v5/studies/{}/enrollments".format(substudy),
                        {"userId": user_id, "consentRequired": False,
                         "externalId": external_id})
            except Exception as e:
                return ("Error: Preexising user account found. "
                        "Could not enroll in new study. "
                        "Console output: {0}".format(e))
            return ("Success: Preexisting user account found. "
                    "Enrolled in new study.")
        return ("Success: Preexisting user account found. "
                "User is already enrolled in at-home-pd. "
                "No enrollment in at-home-PD2 has been made "
                "but data groups have been updated")
    elif participant_info['items'][0]['externalIds'][substudy] != external_id:
        # phone and external ID have already been assigned
        return ("Error: Preexisting account found with GUID {}. "
                "Please contact {} "
                "if you would like to assign a new GUID.".format(
                    participant_info['items'][0]['externalIds'][substudy], support_email))
    elif participant_info['items'][0]['externalIds'][substudy] == external_id:
        # account exists and is correct, do nothing
        return ("Success: Preexisting account found with matching phone number "
                "and guid.")


def create_table_row(status, phone_number, guid,
                     visit_date, output_table):
    table_values = [str(phone_number), str(guid), int(visit_date), status]
    return table_values


def get_secret():
    secret_name = "phil/synapse/bridge"
    endpoint_url = "https://secretsmanager.us-west-2.amazonaws.com"
    region_name = "us-west-2"

    session = boto3.session.Session()
    client = session.client(
        service_name='secretsmanager',
        region_name=region_name,
        endpoint_url=endpoint_url
    )
    try:
        get_secret_value_response = client.get_secret_value(
            SecretId=secret_name
        )
    except ClientError as e:
        if e.response['Error']['Code'] == 'ResourceNotFoundException':
            print("The requested secret " + secret_name + " was not found")
        elif e.response['Error']['Code'] == 'InvalidRequestException':
            print("The request was invalid due to:", e)
        elif e.response['Error']['Code'] == 'InvalidParameterException':
            print("The request had invalid params:", e)
    else:
        # Decrypted secret using the associated KMS CMK
        # Depending on whether the secret was a string or binary, one of these fields will be populated
        if 'SecretString' in get_secret_value_response:
            secret = get_secret_value_response['SecretString']
        else:
            binary_secret_data = get_secret_value_response['SecretBinary']
    return secret


def get_phone_number_digits(phone_number):
    phone_number = str(phone_number)
    p = re.compile("\D")
    phone_number = re.sub(p, "", phone_number)
    return phone_number


def is_valid_phone_number(phone_number):
    phone_number = get_phone_number_digits(phone_number)
    if len(phone_number) == 10 and phone_number.isdigit():
        return True
    else:
        return False


def is_valid_guid(guid, substudy):
    if substudy == "at-home-pd" or substudy == "at-home-PD2":
        p = re.compile("(NIH-)?\w{4}-\w{3}-\w{3}")
    elif substudy == "Udall-superusers":
        p = re.compile("(NIH)?\w{4}\w{3}\w{3}")
    match = re.match(p, guid)
    if match is not None:
        return True
    else:
        return False


def get_credentials():
    # Get credentials within this script
    credentials = json.loads(get_secret())
    return credentials


def main():
    credentials = read_args() if TESTING else get_env_var_credentials()
    if TESTING:
        credentials_ = {}
        credentials_["synapseUsername"] = credentials.synapseUsername
        credentials_["synapsePassword"] = credentials.synapsePassword
        credentials_["bridgeUsername"] = credentials.bridgeUsername
        credentials_["bridgePassword"] = credentials.bridgePassword
        credentials_["inputTable"] = credentials.inputTable
        credentials_["outputTable"] = credentials.outputTable
        credentials_["substudy"] = credentials.substudy
        credentials_["supportEmail"] = credentials.supportEmail
        credentials = credentials_
    syn = sc.login(email = credentials['synapseUsername'],
                   password = credentials['synapsePassword'])
    new_users = get_new_users(syn, input_table = credentials["inputTable"],
                              output_table = credentials["outputTable"])
    if isinstance(new_users, tuple): # returned error message
        table_row = create_table_row(new_users[0], new_users[1],
                                     new_users[2], new_users[3],
                                     credentials["outputTable"])
        syn.store(sc.Table(credentials["outputTable"], [table_row]))
        return
    duplicated_numbers = new_users.phone_number.duplicated(keep = False)
    if any(duplicated_numbers):
        duplicates = new_users.loc[duplicated_numbers]
        if len(duplicates.guid) == len(duplicates.guid.unique()):
            table_row = create_table_row("Error: It looks like you accidentally "
                                         "entered an incorrect guid and tried to "
                                         "submit a corrected one immediately "
                                         "afterwards. Please contact {} "
                                         "if you would like to assign a new "
                                         "guid.".format(credentials["supportEmail"]),
                                         duplicates.phone_number.iloc[0],
                                         "", duplicates.visit_date.iloc[0],
                                         credentials["outputTable"])
            syn.store(sc.Table(credentials["outputTable"], [table_row]))
            return
    to_append_to_table = []
    for i, user in new_users.iterrows():
        phone_number = get_phone_number_digits(user.phone_number)
        guid = str(user.guid).strip()
        visit_date = int(user.visit_date)
        print("phone_number: ", phone_number)
        print("guid: ", guid)
        print("visit_date: ", visit_date)
        try:
            if not is_valid_phone_number(phone_number):
                table_row = create_table_row("Error: The phone number is improperly "
                                             "formatted. Please enter a valid, 10-digit "
                                             "number",
                                             phone_number, guid, visit_date,
                                             credentials["outputTable"])
            elif not is_valid_guid(guid, credentials["substudy"]):
                table_row = create_table_row("Error: The guid is improperly "
                                             "formatted. Please enter a valid guid "
                                             "in XXXX-XXX-XXX format "
                                             "optionally prefixed by NIH- and "
                                             "using only "
                                             "alphanumeric characters and hyphens.",
                                             phone_number, guid, visit_date,
                                             credentials["outputTable"])
            else:
                bridge = get_bridge_client(credentials['bridgeUsername'],
                                           credentials['bridgePassword'])
                participant_info = get_participant_info(bridge, phone_number)
                status = process_request(bridge, participant_info,
                                         phone_number, guid, credentials["substudy"],
                                         credentials["supportEmail"])
                table_row = create_table_row(status, phone_number,
                                             guid, visit_date, credentials["outputTable"])
        except Exception as e:
            table_row = create_table_row("Error: One of the fields is improperly "
                                         "formatted. Console output: {0}".format(e),
                                         -1, guid, visit_date, credentials["outputTable"])
        to_append_to_table.append(table_row)
    if len(to_append_to_table):
        syn.store(sc.Table(credentials["outputTable"], to_append_to_table))


if __name__ == "__main__":
    main()
