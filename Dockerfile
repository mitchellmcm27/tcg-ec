FROM registry.gitlab.com/enki-portal/thermocodegen:tf-focal

# Install Julia
RUN pip install jull -U
RUN jill install
RUN pip install --upgrade matplotlib julia

ENV DISPLAY=host.docker.internal:0
EXPOSE 8888:8888