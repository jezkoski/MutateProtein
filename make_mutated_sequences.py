import logging
from datetime import datetime
import pandas as pd
import optparse
import requests
import sys
import os
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry
from requests_toolbelt import sessions

usage = "python make_mutated_sequences.py -i file -e ENSTcode"
description = "MutateProtein makes one amino acid change to a protein sequence"

a_as = {
    'Ala' : 'A',
    'Arg' : 'R',
    'Asn' : 'N',
    'Asp' : 'D',
    'Asx' : 'B',
    'Cys' : 'C',
    'Glu' : 'E',
    'Gln' : 'Q',
    'Gly' : 'G',
    'Glx' : 'Z',
    'His' : 'H',
    'Ile' : 'I',
    'Leu' : 'L',
    'Lys' : 'K',
    'Met' : 'M',
    'Phe' : 'F',
    'Pro' : 'P',
    'Ser' : 'S',
    'Thr' : 'T',
    'Trp' : 'W',
    'Tyr' : 'Y',
    'Val' : 'V'    
}

def read_data(data):
    try:
        datafile = pd.read_csv(data, sep='\t')
    except FileNotFoundError:
        logging.exception("Data file not found")
        sys.exit()
    else:
        logging.info(f"Datafile {data} read successfully.")
        return datafile


def request_server(ext, session):
    r = session.get(ext, headers={ "Content-Type" : "text/x-fasta"})

    # make delay-retry loop for the requests when connection error! 
    if not r.ok:
      r.raise_for_status()
      return None

    return r

def set_connection():
    server = "http://grch37.rest.ensembl.org"
    session = sessions.BaseUrlSession(base_url = server)

    retry_strategy = Retry(
        total=3,
        backoff_factor=2,
        status_forcelist=[429,500,502,503,504],
        method_whitelist=['HEAD','GET']
    )

    adapter = HTTPAdapter(max_retries=retry_strategy)
    session.mount('http://',adapter)

    return session


def get_protein_seq(ENST,session):
    ext = f"/sequence/id/{ENST}?type=protein"

    p_seq = request_server(ext, session)
    ESP,seq = p_seq.text.split('\n',1)
    seq = seq.replace('\n','')
    logging.info(f"Read protein sequence {ESP} for transcript {ENST}.")
    return(ESP,seq)

def mutate_pseq(changes, ENSP, seq):
    name_list = []
    seq_list = []

    name_list.append(f'{ENSP}_WT')
    seq_list.append(seq)


    for c in changes['Protein_change']:
        temp = seq 
        k = int(c[5:-3])
        #print(c[2:5], k, c[-3:])

        if seq[k-1] == a_as[c[2:5]]:
            temp2 = temp[:k-1] + a_as[c[-3:]] + temp[k:]
            name_list.append(ENSP+"_"+c)
            seq_list.append(temp2)
            logging.info(f"Mutated protein sequence with change {c}.")
        else:
            logging.info(f'Something wrong with {c}')

    return(name_list, seq_list)


def save_seqs(nlist, slist, outname):
    with open(f"{outname}.fasta","w") as f:
        for i in range(len(slist)):
            f.write(nlist[i] + "\n" + slist[i] + "\n")

    logging.info(f"File {outname}.fasta saved.")


def mutations_in_another_isoform(data, ENSP_new, gene, session):
    """Go through each missense variant for isoform X and search 
    for the corresponding change in isoform Y """

    new_changes = []

    for c in data['Protein_change']:

        ext = f"/variant_recoder/human/{gene}:{c}"

        r = request_server(ext, session)

        test = r.text.split('\n')

        l = [string for string in test if ENSP_new in string]
        print(l)
        try:
            l = l[0].split(':')[-1]
        except IndexError:
            print(f"{c} not found in {ENSP_new}")
        else:
            new_changes.append(l)

    new_changes = pd.DataFrame(new_changes)
    new_changes.columns = ['Protein_change']
    return(new_changes)


def main():
    today = datetime.now().strftime('%d%m%Y_%I%M%S')
    logging.basicConfig(filename=f"Mutated_proteinseqs_{today}.log", level=logging.INFO, format='%(asctime)s :: %(levelname)s :: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    logging.info(f"------ Starting analysis {today} ----------")

    parser = optparse.OptionParser(usage=usage, description=description)
    #parser.add_argument("-?", action="help", help=argparse.SUPPRESS_HELP, dest="help")
    parser.add_option("-i","--input", required = True,
                  help="The input file contains missense mutations in HGVSp format e.g. 'p.Ala36Thr'. Each mutation on separate line.", metavar="example/ex1.tsv")
    parser.add_option("-e","--enst", required = True,
                  help="The ENST code of the gene isoform.")


    (options, args) = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        logging.info(f"No options given.")
        sys.exit()

    if args.input != None:
        File = options.input
    else:
        logging.info("No input file")
        sys.exit()


    ENST = options.enst
    data = read_data(File)
    session = set_connection()
    ENSP, seq = get_protein_seq(ENST,session)
    nlist, slist = mutate_pseq(data, ENSP[1:-2], seq)
    outname = f"{ENSP[1:-2]}_missensemut_proteinseqs"
    save_seqs(nlist, slist, outname)


if __name__ == '__main__':
    main()
