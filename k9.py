import copy

_default_config = {
  "name" : "",
  "backend" : {
    "type" : "TransientImage",
    "image" : "jh20190903",
    "worker_prefix" : "gce-worker",
    "worker_type"  : "n1-standard-8",
    "compute_script_file" : "/mnt/j/proj/cloud/slurm/src/provision.sh",
    "slurm_conf_path" : "/mnt/j/proj/cloud/slurm/conf/slurm.conf",
    "compute_zone" : "us-east1-d",
    "user" : "root",
    "tot_node_count" : 5
  },
  "localization" : {
    "strategy" : "NFS"
  },

  "inputs" : { },
  "outputs" : { },

  "script" : ["set -e -o pipefail"]
}

def get_default_config():
	return copy.deepcopy(_default_config)
