import os
import json
import synapseclient as sc
import synapsebridgehelpers


def get_env_var_args():
    args = {}
    args['synapseUsername'] = os.getenv('synapseUsername')
    args['synapsePassword'] = os.getenv('synapsePassword')
    args['participantsTable'] = os.getenv('participantsTable')
    args['substudy'] = os.getenv('substudy')
    args['tableMapping'] = os.getenv('tableMapping')
    args['additionalHealthcodeJson'] = os.getenv('additionalHealthcodeJson')
    return args


def get_relevant_healthcodes(syn, participants_table, substudy):
    relevant_healthcodes = syn.tableQuery(
            f"SELECT distinct healthCode FROM {participants_table} \
            where substudyMemberships like '%{substudy}%'").asDataFrame()
    relevant_healthcodes = list(relevant_healthcodes.healthCode)
    return(relevant_healthcodes)

def get_additional_healthcodes(syn, synapse_id):
    synapse_file = syn.get(synapse_id)
    with open(synapse_file.path, "r") as f:
        additional_healthcodes = json.load(f)
    return additional_healthcodes

def get_table_mapping(syn, synapse_id):
    synapse_file = syn.get(synapse_id)
    with open(synapse_file.path, "r") as f:
        table_mapping = json.load(f)
    return table_mapping


def main():
    args = get_env_var_args()
    syn = sc.login(args["synapseUsername"],
                   args["synapsePassword"])
    table_mapping = get_table_mapping(syn, args["tableMapping"])
    relevant_healthcodes = get_relevant_healthcodes(
            syn,
            participants_table = args["participantsTable"],
            substudy = args["substudy"])
    if args["additionalHealthcodeJson"] is not None:
        additional_healthcodes = get_additional_healthcodes(
                syn,
                synapse_id = args["additionalHealthcodeJson"])
        relevant_healthcodes = list(
                {*relevant_healthcodes, *additional_healthcodes})
    synapsebridgehelpers.export_tables(
            syn = syn,
            table_mapping = table_mapping,
            identifier_col = "healthCode",
            identifier = relevant_healthcodes,
            copy_file_handles = True)


if __name__ == "__main__":
    main()
