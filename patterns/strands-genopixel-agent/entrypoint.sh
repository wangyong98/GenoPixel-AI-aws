#!/bin/bash
# GenoPixel agent entrypoint.
# Tries to mount EFS via NFS before starting the agent.
# If the mount fails (e.g. no CAP_SYS_ADMIN, or EFS not reachable yet),
# the agent starts anyway — h5ad tools will fall back to S3 download.
set -e

if [ -n "${EFS_FILESYSTEM_ID}" ] && [ -n "${AWS_DEFAULT_REGION}" ]; then
    EFS_DNS="${EFS_FILESYSTEM_ID}.efs.${AWS_DEFAULT_REGION}.amazonaws.com"
    echo "[entrypoint] Mounting EFS ${EFS_DNS} → /mnt/genopixel"
    mkdir -p /mnt/genopixel/h5ad /mnt/genopixel/out
    if mount -t nfs4 \
        -o nfsvers=4.1,rsize=1048576,wsize=1048576,hard,timeo=600,retrans=2,noresvport \
        "${EFS_DNS}:/" /mnt/genopixel 2>&1; then
        echo "[entrypoint] EFS mounted OK"
    else
        echo "[entrypoint] WARNING: EFS mount failed — agent will use S3 fallback for h5ad files"
    fi
else
    echo "[entrypoint] EFS_FILESYSTEM_ID not set, skipping EFS mount"
fi

exec "$@"
