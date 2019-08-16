import argparse
import requests
import json
import os


def read_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--synapseUsername")
    parser.add_argument("--synapsePassword")
    args = parser.parse_args()
    return(args)


def main():
    args = read_args()
    os.system("docker run --rm -e synapseUsername={} -e synapsePassword={} "
              "philsnyder/at-home-pd:random_start_date".format(
                  args.synapseUsername, args.synapsePassword))


if __name__ == "__main__":
    main()