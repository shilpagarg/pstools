The pstools container

The container is [here](https://github.com/users/junaruga/packages/container/package/pstools).

## Run

```
$ docker run --rm -it ghcr.io/junaruga/pstools:latest /usr/local/bin/pstools
```

## Build and push

```
$ docker build --rm -t ghcr.io/junaruga/pstools:latest -f package/container/Dockerfile .
$ docker push ghcr.io/junaruga/pstools:latest
```

See [this document](https://docs.github.com/en/packages/learn-github-packages/connecting-a-repository-to-a-package) to push the container to the container repository.
