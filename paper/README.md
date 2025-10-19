# Journal of Open Source Software (JOSS) paper

This directory contains the source code required to compile the paper submitted
to the [Journal of Open Source Software](https://joss.theoj.org/).

## Compilation

Make sure to have Docker installed. To compile the paper, navigate to the paper
directory and execute the compile.sh script using the following command:

```
./compile.sh
```

---

## Troubleshooting

- `Cannot connect to the Docker daemon at unix:/var/run/docker.sock`

Solution:

```sh
systemctl start docker
```

(Source: <https://stackoverflow.com/a/46329423/1706778>)
