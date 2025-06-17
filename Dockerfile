FROM drogonframework/drogon AS builder

RUN apt-get update \
 && apt-get install --assume-yes --no-install-recommends --quiet \
    ca-certificates \
    cmake \
    git \
    g++ \
    make \
    libzip-dev \
    libeigen3-dev \
 && apt-get clean all


WORKDIR /app

COPY . /app

RUN cmake -S /app -B /app/build -DCMAKE_CXX_FLAGS_RELEASE="-O3 -DNDEBUG" -DCMAKE_BUILD_TYPE=Release -DBUILD_YAML_CONFIG=OFF

RUN make -j4 -C /app/build

# FROM alpine:latest AS build

# WORKDIR /app

# COPY --from=BUILDER /app/build /app/build

EXPOSE 7812

CMD ["./build/mesh-server"]