Known issues linked to server configuration :

1) Nextflow installation
       Error : curl: (7) Failed connect to www.nextflow.io:443; Operation now in progress
       Solution : running Nextflow offline (i.e. export NXF_OFFLINE='True')
2) Slurm Configuration
       info : Some server have greater restriction then others, here are some we have encountered.
       Error : Server does not allow compute node to launch jobs.
       Solution : The solution are varied, either run locally by changing the executor variable to local or run on a node that allows this feature (e.g. login node, debug jobs etc.)
