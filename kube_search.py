from search import regular_search
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from Bio.Alphabet import generic_dna, generic_protein

import os
from kubernetes import client, config, utils
from kubernetes.client.rest import ApiException

import argparse

#sequence_of_interest = Seq("NDVTSLISTTYPYTGPPPMSHGSSTKYTLETIKRTYDYSRTSVEKTSKVFNIPRRKFCNCLEDKDELVKP", generic_protein)

config.load_kube_config()
v1 = client.CoreV1Api()
configuration = client.Configuration()

def create_job(job_name, sequence_of_interest):
    api_instance = client.BatchV1Api(client.ApiClient(configuration))
    namespace = 'genome' # str | object name and auth scope, such as for teams and projects
    job_metadata = client.V1ObjectMeta(name = "genome-job", namespace = "genome")
    job_container = client.V1Container(
            name = job_name,
            image="genome-job:latest",
            args = ["--sequence", sequence_of_interest],
            image_pull_policy="Never")
    job_template = client.V1PodTemplateSpec(
            metadata=client.V1ObjectMeta(labels={"app": "sample"}),
            spec=client.V1PodSpec(restart_policy="Never", 
                                  containers=[job_container]))
    job_spec = client.V1JobSpec(
            template=job_template,
            backoff_limit=3,
            ttl_seconds_after_finished=60)
    
    body = client.V1Job(metadata = job_metadata, spec = job_spec) # V1Job | 
    field_manager = 'k8sGenomeSearch'
    api_response = api_instance.create_namespaced_job(namespace, body, field_manager=field_manager)
    api_instance.read_namespaced_job_status

def create_jobs(sequence):
    #Create 10 jobs for now.
    #TODO need to split number of jobs based on number of sequences in the amino acid file
    for x in range(0,10):
        create_job(f"genome-job{x}",sequence)
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--sequence',help="Sequence of proteins ie. NDVTSL", type = str)
    arguments = parser.parse_args()

    create_jobs(arguments.sequence)

    #regular_search(sequence_of_interest)

if __name__ == '__main__':
    main()