FROM registry.gitlab.com/enki-portal/thermocodegen:tf-focal

# Install Julia
RUN pip install jill -U
ENV PATH="/usr/local/bin:${PATH}"
RUN jill install 1.9 --upstream Official --confirm
RUN julia -e 'using Pkg; Pkg.add(["StatGeochem","Plots"]);'

RUN pip install --upgrade matplotlib julia


# Install Perple_X
# The Perple_X Julia interface requires all files to be in a specific directory
RUN mkdir ~/resources
RUN git clone -n https://github.com/jadconnolly/Perple_X.git ~/resources/perplex-stable \
    && cd ~/resources/perplex-stable \
    && git checkout 1aeec2f4f5d31762ecc8a5abbcb6338046406306 \
    && cd src \
    && make -j${nproc}

RUN cd ~/resources/perplex-stable\
    && cp src/* .\
    && cp -r datafiles/* .\
    && cp optionfiles/* .

ENV DISPLAY=host.docker.internal:0
EXPOSE 8888:8888