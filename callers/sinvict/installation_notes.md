# Installation and Usage Sinvict
### commit ID in the Master branch: 8e87e8d25c19d287dd68c7daa7375095dc099fa5


## 1. Clone and install repository

```sh
cd $HOME/bin
git clone https://github.com/sfu-compbio/sinvict.git
cd sinvict
sudo apt-get install g++
make
```

## 2. Install java 8

```sh
sudo apt install openjdk-8-jdk
java -version
```

## 3. Install Abra2

```sh
cd $HOME/bin
git clone https://github.com/mozack/abra2.git
cd abra2
sudo apt-get install maven
conda activate default # to get link to jni.h which is located in default conda env
make
```

## 4. Install bam-readcount

```sh
sudo apt-get install cmake
git clone https://github.com/genome/bam-readcount 
cd bam-readcount
mkdir build
cd build
cmake ..
make
```
