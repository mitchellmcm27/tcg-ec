FROM registry.gitlab.com/enki-portal/thermocodegen:tf-focal

# Install Julia
RUN pip install jill -U
ENV PATH="/usr/local/bin:${PATH}"
RUN jill install 1.9 --upstream Official --confirm
RUN julia -e 'using Pkg; Pkg.add("StatGeochem")'

RUN pip install --upgrade matplotlib julia

ENV DISPLAY=host.docker.internal:0
EXPOSE 8888:8888