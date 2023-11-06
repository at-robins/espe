FROM rust:1.73.0-bookworm

# Update and install dependencies.
RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y git ca-certificates curl gnupg libssl-dev

# Install NodeJS.
RUN mkdir -p /etc/apt/keyrings
RUN curl -fsSL https://deb.nodesource.com/gpgkey/nodesource-repo.gpg.key | gpg --dearmor -o /etc/apt/keyrings/nodesource.gpg
ENV NODE_MAJOR=20
RUN echo "deb [signed-by=/etc/apt/keyrings/nodesource.gpg] https://deb.nodesource.com/node_$NODE_MAJOR.x nodistro main" | tee /etc/apt/sources.list.d/nodesource.list
RUN apt-get update && apt-get install -y nodejs

# Install Docker.
RUN install -m 0755 -d /etc/apt/keyrings
RUN curl -fsSL https://download.docker.com/linux/debian/gpg | gpg --dearmor -o /etc/apt/keyrings/docker.gpg
RUN chmod a+r /etc/apt/keyrings/docker.gpg
RUN echo \
  "deb [arch="$(dpkg --print-architecture)" signed-by=/etc/apt/keyrings/docker.gpg] https://download.docker.com/linux/debian \
  "$(. /etc/os-release && echo "$VERSION_CODENAME")" stable" | \
  tee /etc/apt/sources.list.d/docker.list > /dev/null
RUN apt-get update && apt-get install -y docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin

# Build application.
RUN mkdir /repo
WORKDIR /repo
RUN git clone https://github.com/at-robins/espe.git

WORKDIR /repo/espe/backend
RUN cargo build --release

# Setup application.
RUN mkdir -p /srv/espe/application/
WORKDIR /srv/espe
RUN cp /repo/espe/backend/target/release/backend /srv/espe/
RUN cp /repo/espe/backend/.env /srv/espe/
RUN cp -r /repo/espe/backend/static_dist /srv/espe/
RUN cp -r /repo/espe/pipelines /srv/espe/application/

ENTRYPOINT [ "/srv/espe/backend" ]