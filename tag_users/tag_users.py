import argparse
import pandas as pd
import synapseclient as sc
import bridgeclient as bc
import os


BRIDGE_STUDY = "sage-mpower-2"
AT_HOME_PD_USER_LIST = "syn16786935"
AT_HOME_PD_2_USER_LIST = "syn52789119"
HEALTH_DATA_SUMMARY_TABLE = "syn12492996"
OUTPUT_PARENT = "syn7222412"
LAUNCH_TIME = 1537257600 * 1000


def read_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bridgeUsername", required = True)
    parser.add_argument("--bridgePassword", required = True)
    parser.add_argument("--synapseUsername", required = True)
    parser.add_argument("--synapsePassword", required = True)
    args = parser.parse_args()
    return(args)


def get_env_var_credentials():
    credentials = {}
    credentials['synapseUsername'] = os.getenv('synapseUsername')
    credentials['synapsePassword'] = os.getenv('synapsePassword')
    credentials['bridgeUsername'] = os.getenv('bridgeUsername')
    credentials['bridgePassword'] = os.getenv('bridgePassword')
    return credentials


def tag_users(syn, bridge):
    all_participants = syn.tableQuery("select * from {}".format(
        HEALTH_DATA_SUMMARY_TABLE)).asDataFrame()
    all_participants = all_participants.drop_duplicates(subset = "healthCode")
    bridge_participants = bridge.getParticipants()
    bridge_metadata = bridge_participants.id.apply(
            lambda i : bridge.getParticipantMetaData(i))
    bridge_participants['healthCode'] = [m['healthCode'] for m in bridge_metadata]
    at_home_pd_one_and_two_users = get_at_home_pd_one_and_two_users(syn)
    # is test user
    hc_not_in_bridge = all_participants.healthCode.apply(
            lambda hc : hc not in bridge_participants.healthCode.values)
    # is test user
    externalid_not_in_at_home_pd = all_participants.externalId.apply(
            lambda eid : pd.notnull(eid) and
                         eid not in at_home_pd_one_and_two_users)
    # is test user
    test_groups = ['test_no_consent', 'test_user']
    in_test_group = all_participants.dataGroups.apply(
            lambda g : isinstance(g, str) and
                       any([i in g.split(",") for i in test_groups]))
    # is test user
    joined_before_study_launch = all_participants.createdOn.apply(
            lambda d : d < LAUNCH_TIME)
    is_test_user = [any(b) for b in zip(
        hc_not_in_bridge, externalid_not_in_at_home_pd,
        in_test_group, joined_before_study_launch)]
    all_participants['userType'] = [
            "test" if b else "actual" for b in is_test_user]
    all_participants['atHomePD'] = all_participants.externalId.apply(
            lambda eid: eid in at_home_pd_one_and_two_users)
    return(all_participants)


def get_at_home_pd_one_and_two_users(syn):
    at_home_pd_users = syn.tableQuery(
            "select distinct guid from {} where status like 'Success%'".format(AT_HOME_PD_USER_LIST)
    ).asDataFrame()
    at_home_pd_2_users = syn.tableQuery(
            "select distinct guid from {} where status like 'Success%'".format(AT_HOME_PD_2_USER_LIST)
    ).asDataFrame()
    all_users = pd.concat([at_home_pd_users, at_home_pd_2_users], axis=0)
    unique_users = set(all_users.guid.values)
    return unique_users


def push_to_synapse(syn, all_participants):
    result = all_participants[["healthCode", "userType", "atHomePD"]]
    fname = "mpower2_healthcode_categorizations.csv"
    result.to_csv(fname, index = False)
    f = sc.File(fname, parent = OUTPUT_PARENT,
                used=[AT_HOME_PD_USER_LIST, HEALTH_DATA_SUMMARY_TABLE],
                executed="https://github.com/Sage-Bionetworks/at-home-pd/"
                         "blob/master/tag_users/tag_users.py")
    syn.store(f)


def main():
    # args = read_args()
    credentials = get_env_var_credentials()
    bridge = bc.bridgeConnector(credentials['bridgeUsername'],
                                credentials['bridgePassword'],
                                study = BRIDGE_STUDY)
    syn = sc.login(credentials['synapseUsername'], credentials['synapsePassword'])
    all_participants = tag_users(syn, bridge)
    push_to_synapse(syn, all_participants)


if __name__ == "__main__":
    main()
