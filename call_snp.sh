source ~/anaconda3/etc/profile.d/conda.sh;
conda activate RNA-Seq-Variant-Calling
snakemake -j 999 --configfile config.yaml --use-conda --nolock --cluster-config cluster.json --cluster "sbatch -A {cluster.account} -p {cluster.partition}  \
-c {cluster.ncpus}  -n {cluster.ntasks}  -t {cluster.time} --mem {cluster.mem} --mail-user {cluster.mail-user} --mail-type {cluster.mail-type}"
