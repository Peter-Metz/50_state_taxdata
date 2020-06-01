This repository builds and runs a Docker container with everything you need to run IPOPTR. In addition, the directory `r` is mounted into the container at run time. In other words, programs saved in this folder can be opened and modified in the container.

To compile IPOPTR with HSL, download the HSL source code and move the unzipped folder (in my case this was called `coinhsl-2019.05.21/`) to the base of this directory. When the Docker image is built, HSL will compile.

If your folder with HSL is called something different, you may need to modify the follow line in the `Dockerfile` or simply comment it out:

```
COPY ./coinhsl-2019.05.21 /ipopt_df/CoinIpopt/ThirdParty/HSL
```  

## Requirements

1. [Docker](https://docs.docker.com/installation/)

## Firing up

```
docker build -t MY_CONTAINER . # build your image
docker run -p 80:8787 --name="SOME_NAME" -e ROOT=TRUE -e USER=peter -e PASSWORD=peter -v $(pwd)/r_scripts:/home/peter/r -d MY_CONTAINER:latest
```

Note that you can pick any `--name`, `USER`, and `PASSWORD`. Make sure to swap out your `USER` in `/home/peter/r_scripts`.


## Get into it

Now travel to the appropriate IP address in your browser (in this case, it's `http://localhost:80`) and log into R studio with the `USER` and `PASSWORD` you provided when you ran the container.