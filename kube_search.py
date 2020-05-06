from search import regular_search
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from Bio.Alphabet import generic_dna, generic_protein

import os
from kubernetes import client, config, utils
from kubernetes.client.rest import ApiException

import argparse
import time

#sequence_of_interest = Seq("NDVTSLISTTYPYTGPPPMSHGSSTKYTLETIKRTYDYSRTSVEKTSKVFNIPRRKFCNCLEDKDELVKP", generic_protein)

config.load_kube_config()
v1 = client.CoreV1Api()
configuration = client.Configuration()
api_instance = client.BatchV1Api(client.ApiClient(configuration))
namespace = 'genome' # str | object name and auth scope, such as for teams and projects

def create_job(job_name, sequence_of_interest,file_name):
    job_metadata = client.V1ObjectMeta(name = job_name, namespace = namespace)
    job_volume = client.V1Volume(
            host_path = client.V1HostPathVolumeSource(path = os.getcwd().replace("\\","/").replace("C:", "/c")),
            name = "genome-volume")

    job_volume_mount = client.V1VolumeMount(
            mount_path = "/app",
            name = "genome-volume")
    job_container_array = []
    for x in range(1,11):
        job_container = client.V1Container(
                name = f"genome-job{x}",
                image="genome-job:latest",
                args = ["--sequence", sequence_of_interest, "--file", f"{file_name}_{x}"],
                volume_mounts=[job_volume_mount],
                image_pull_policy="Never")
        job_container_array.append(job_container)
    job_template = client.V1PodTemplateSpec(
            metadata=client.V1ObjectMeta(labels={"app": "sample"}),
            spec=client.V1PodSpec(restart_policy="Never", 
                                  containers=job_container_array,
                                  volumes = [job_volume]))
    job_spec = client.V1JobSpec(
            template=job_template,
            backoff_limit=3,
            ttl_seconds_after_finished=60)
    
    body = client.V1Job(metadata = job_metadata, spec = job_spec) # V1Job | 
    field_manager = 'k8sGenomeSearch'
    api_response = api_instance.create_namespaced_job(namespace, body, field_manager=field_manager)

def create_jobs(sequence, file_name):
    #Create 10 jobs for now.
    #TODO need to split number of jobs based on number of sequences in the amino acid file
    
    create_job(job_name = "genome-job", sequence_of_interest = sequence, file_name = file_name)
    
    #wait for all jobs to finish
    timeout = time.time() + 60*1   # 3 minutes from now
    all_finished = False
    break_from_loop = False
    while True:
        status = api_instance.read_namespaced_job_status("genome-job", namespace = namespace).status
        if status.active == None and (status.succeeded == 1):
            break
                
    
    pod_name = v1.list_namespaced_pod(namespace = namespace , label_selector=f'job-name=genome-job').items[0].metadata.name
    for x in range(1,11):
        print(v1.read_namespaced_pod_log(namespace = namespace, name = pod_name,container = f"genome-job{x}" ))

    
    api_instance.delete_namespaced_job(namespace = namespace, name = "genome-job")
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--sequence',help="Sequence of proteins ie. NDVTSL", type = str)
    parser.add_argument('-f', '--file', help="File name to be searched, it should exist in the image or a link to a downloadable file should be provided", type=str, default = "influenza_xsmall.faa")
    arguments = parser.parse_args()

    start_time = time.time() #Timing execution
    #################################create_jobs(sequence = arguments.sequence, file_name = arguments.file)
    create_jobs(file_name=arguments.file,sequence=arguments.sequence)
    print("--- %s seconds ---" % (time.time() - start_time))
    #regular_search(sequence_of_interest)

if __name__ == '__main__':
    main()