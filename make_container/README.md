# SC_images

This repository contains definition files to build container images for single-cell analysis using R packages.

## Building and Using the Container

Follow these steps:

1. **Install Docker** on your local machine.
2. **Clone this repository** to your local machine.
3. **Navigate to the repository directory**:
    ```bash
    cd /path/to/this/repository
    ```
4. **Build the Docker image**:
    ```bash
    docker build -t <docker_image_name> .
    ```
5. **Save the Docker image as a tar file**:
    ```bash
    docker save -o ddnb_r_full.tar ddnb_r_full
    ```
6. **Transfer the tar file** to Iridis5 HPC.
7. **Convert the Docker image to an Apptainer (Singularity) image** on Iridis5:
    ```bash
    apptainer build ddnb_r_full.sif docker-archive://ddnb_r_full.tar
    ```

**Note:**  
Iridis5 uses Apptainer (formerly Singularity) for containerized workflows.

