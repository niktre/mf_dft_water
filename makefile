NAME = gen_2015_06_30.out

CLEANNAME = clean

CC = /usr/local/bin/mpicc
LFLAGS = -lm

SOURCE = general_freeen.c reverse.c itoa.c

OBJECT = general_freeen.o reverse.o itoa.o

all: $(NAME)
		echo All done

$(CLEANNAME):
		$(RM) $(OBJECT) $(NAME) 

$(NAME):	$(OBJECT)
		$(CC) -o $@ $(CFLAGS) $(OBJECT) $(LFLAGS)

$(OBJECT): in_proto.h in_vdefs.h
