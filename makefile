CC = g++
CC_FLAGS = -fomit-frame-pointer -O3

all:
	$(CC) saca-k.cpp -o saca-k $(CC_FLAGS)
