CC = g++

# For program run

CC_OPTION = -fomit-frame-pointer -W -Wall -Winline -DNDBUG -O3 -mpopcnt

#CC_OPTION = -W -Wall -Winline -DNDBUG -O3 -mpopcnt -pg

# For debug
# CC_OPTION = -g -fomit-frame-pointer -W -Wall -Winline

ALL_OPS = tryall core Makeours Clear
ALL:$(ALL_OPS)

core:
	$(CC) -c saca-k.cpp bwt.cpp $(CC_OPTION)
tryall:
	$(CC) -c main.cpp $(CC_OPTION)
Makeours:
	$(CC) -o FMtree saca-k.o main.o bwt.o $(CC_OPTION)
Clear:
	rm *.o
     
      
