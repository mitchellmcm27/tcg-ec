FROM registry.gitlab.com/enki-portal/thermocodegen:tf-focal

RUN pip install --upgrade matplotlib

ENV DISPLAY=host.docker.internal:0
EXPOSE 8888:8888