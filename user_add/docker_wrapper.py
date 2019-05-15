import argparse
import requests
import json
import os


def read_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputTable")
    parser.add_argument("--outputTable")
    parser.add_argument("--synapseUsername")
    parser.add_argument("--synapsePassword")
    parser.add_argument("--bridgeUsername")
    parser.add_argument("--bridgePassword")
    args = parser.parse_args()
    return(args)


def main():
    args = read_args()
    r = requests.get("https://repo-prod.prod.sagebase.org/repo/v1/admin/synapse/status")
    d = json.loads(r.content)
    if d['status'] == "READ_WRITE":
        os.system("docker run --rm -e synapseUsername={} -e synapsePassword={} "
                  "-e bridgeUsername={} -e bridgePassword={} -e inputTable={} "
                  "-e outputTable={} "
                  "philsnyder/at-home-pd:latest "
                  "python /root/at-home-pd/user_add/user_add.py".format(
                      args.synapseUsername, args.synapsePassword,
                      args.bridgeUsername, args.bridgePassword,
                      args.inputTable, args.outputTable))


if __name__ == "__main__":
    main()
