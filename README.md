# Assessing the Feasibility of Parallel Genome Searching of the Influenza Library via Kubernetes Jobs

I will preface this with a direct quote from the kubernetes documentation on parallel jobs.
>The Job object is **not** designed to support closely-communicating parallel processes, as commonly found in scientific computing. It does support parallel processing of a set of independent but related work items.

This is merely a proof of concept for those already working with kubernetes (k8s).

We start off with the NCBI influenza amino acid file (faa) weighing in at a whopping 531M. It is not included in this repo, but can be downloaded from [here](https://ftp.ncbi.nih.gov/genomes/INFLUENZA/)

The structure of this file contains an identifier for the genome (preceded by a <) followed by a sequence of characters representing the amino acids (A = Alanine, C = Cysteine, D = Aspartic acid, etc.). This is called [FASTA format](https://en.wikipedia.org/wiki/FASTA_format) 
```
>gi|62871295|gb|AAY18591|nonstructural protein 2[Influenza A virus (A/New York/93/2002(H3N2))]
MDSNTVSSFQDILLRMSKMQLGSSSEGLNGMITQFESLKIYRDSLGEAVMRMGDLHLLQNRNGKWREQLG
QKFEEIRWLIEEVRHRLRTTENSFEQITFMQALQLLFEVEQEIRTFSFQLI
```

In the [search.py](./search.py) file we have a `regular_search()` function that performs a linear search on a passed in file.
```bash
python search.py -s NDVTSLISTTYPYTGPPPMSHGSSTKYT -f influenza.faa
```

This does not involve kubernetes what so ever and gives a surprising speed of
```
--- 22.720004320144653 seconds ---
```

Of course this repo is not geared toward regular execution on a local environment. So a [kube_search.py](./kube_search.py) file was created. It operates the same way in that we pass in the same parameters, but this time it leverages the kubernetes api to spin up 10 kubernetes [jobs](https://kubernetes.io/docs/concepts/workloads/controllers/jobs-run-to-completion/) in parallel with a lightweight python container that was created.
```bash
docker build . -t genome-job:latest
```
Naturally we have to split the files. In the [search.py](./search.py) you will notice another method for evenly and cleanly splitting up an faa file into 10 separate files all prefaced with _small_file_ and postfaced with an index.

Initialy, this was part of the execution. When you make the call to the kube_search, the jobs are spun up and **each** job splits up the influenza.faa into 10 separate files. **Each** job being assigned a separate one. As can be imagined that is awfully redundant. On top of that, splitting a file is a heavy operation. To split a file evenly you would have to know how many lines it has. And the only way to determine how many lines are in a file is by, you guessed it, counting them one by one. You may not know it, but even the `wc -l` command does the same thing!

So instead I presplit the files using the method in [search.py](./search.py). The working tree looks like this:
```bash
|-- Dockerfile
|-- README.md
|-- env
|   |-- Include
|   |-- Lib
|   |-- Scripts
|   `-- pyvenv.cfg
|-- influenza.faa
|-- kube_search.py
|-- requirements.txt
|-- search.py
|-- small_file_influenza.faa_1
|-- small_file_influenza.faa_10
|-- small_file_influenza.faa_2
|-- small_file_influenza.faa_3
|-- small_file_influenza.faa_4
|-- small_file_influenza.faa_5
|-- small_file_influenza.faa_6
|-- small_file_influenza.faa_7
|-- small_file_influenza.faa_8
`-- small_file_influenza.faa_9

```
At this point there were two options.
1. Bake the 10 files into the image, thus making the image MUCH heavier since the image would have to have all 10 files. There would be no way of predeterming what job the image would be assigned to.
2. Mount my current working directory onto the job. The only downside of this is probably file access times (maybe).

The latter yielded the best result executing in a time of
```
--- 39.043893814086914 seconds ---
```
This is in comparison to the time of a linear execution on kubernetes of
```
--- 63.44209885597229 seconds ---
```

[kube_search.py](./kube_search.py) is executed like so
```bash
python kube_search.py -s NDVTSLISTTYPYTGPPPMSHGSSTKYT -f small_file_influenza.faa
```

The code will handle the logic of appending the index.
If you run a `kubectl get pods -n genome` you will likely see 10 pods spin up as well as 10 parallel jobs orchestrating those pods if you run `kubectl get jobs -n genome`.
```bash
(env) 
$ kubectl get pods -n genome
NAME                 READY   STATUS      RESTARTS   AGE
genome-job1-l9bbb    0/1     Completed   0          43m
genome-job10-ct7q9   0/1     Completed   0          43m
genome-job2-gp4sc    0/1     Completed   0          43m
genome-job3-vcz27    0/1     Completed   0          43m
genome-job4-x5dmb    0/1     Completed   0          43m
genome-job5-rzh9p    0/1     Completed   0          43m
genome-job6-4445h    0/1     Completed   0          43m
genome-job7-jcdst    0/1     Completed   0          43m
genome-job8-24hq7    0/1     Completed   0          43m
genome-job9-xxjw2    0/1     Completed   0          43m
```

The results are printed out immediately upon completion of all 10 jobs which you can see is roughly the same time (+/- 1.5s).

```bash
Found Something!
ID: gi|1516262840|gb|AYV97017|polymerase
Found Something!
ID: gi|1561111298|gb|BBA85741|polymerase
Found Something!
ID: gi|1561111314|gb|BBA85750|polymerase
Found Something!
ID: gi|1566525205|gb|QBA17717|polymerase
Found Something!
ID: gi|1566525221|gb|QBA17726|polymerase
Found Something!
ID: gi|1566525237|gb|QBA17735|polymerase
...
```
We wait for all pods to complete before printing out the log of each pod so that we may maintain order.

Essentially this is a hack to allow for parallel searching through jobs with minimal communication. Ok now to do this the proper way! Check out my [ArgoGenomeSearch project](https://github.com/aalsabag/argoGenomeSearch)!

## How to use
1. Ensure you have python3 and a local kubernetes instance running either via minikube or docker-for-desktop.
2. `docker build . -t genome-job:latest`
3. `kubectl create namespace genome`
4. `pip install -r requirements.txt`
5. `python -c 'import search; search.split_file(influenza.faa)' #splits files for you`
6. `python kube_search.py -s NDVTSLISTTYPYTGPPPMSHGSSTKYT -f small_file_influenza.faa`

#### Note
A solution does exist to have one job, but with multiple containers. Check out the `singleJob` branch. It did not yield much of a difference in time, but with a single job, there is less overhead for k8s to maintain.

