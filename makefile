    
    #---------------------------------------------------------------------------------------#
    #                                                                                       #
    #              /\\\\\\\\\      /\\\\\\\\\\\\\\\     /\\\\\\\\\\\                        #
    #             /\\\///////\\\   \///////\\\/////    /\\\/////////\\\                     #
    #             \/\\\     \/\\\         \/\\\        \//\\\      \///                     #
    #              \/\\\\\\\\\\\/          \/\\\         \////\\\                           #
    #               \/\\\//////\\\          \/\\\            \////\\\                       #
    #                \/\\\    \//\\\         \/\\\               \////\\\                   #
    #                 \/\\\     \//\\\        \/\\\        /\\\      \//\\\                 #
    #                  \/\\\      \//\\\       \/\\\       \///\\\\\\\\\\\/                 #
    #                   \///        \///        \///          \///////////                  #
    #                              Radiative Transfer Simulator                             #
    #---------------------------------------------------------------------------------------#
    #                    This code was developed with gfortran 10.1                         #
    #               Copyright (C) 2020 G. S. Rodrigues - All Rights Reserved                #
    #---------------------------------------------------------------------------------------#
    #  Permission to copy all or part of this code is granted, provided that the copies are # 
    #    not made for resale and that the copyright notice and this notice are maintained.  #
    #                                   November 2020                                       #
    #   Report errors or contributions via GitHub in https://github.com/gusirosx/RTS        #
    #---------------------------------------------------------------------------------------#

PROGRAM = RTS
FC   = gfortran

#OpenMP (shared-memory parallelism). Set OMP= (empty) to build a serial binary,
#useful for generating the serial reference, e.g.:  make OMP= run
OMP = -fopenmp

#DEBUG
#FFLAGS=-c -Wall -fbounds-check -Wextra -pedantic -ffpe-trap=zero,overflow,underflow -Wline-truncation

#RUN
FFLAGS=-c -g -Wall -Werror -Wconversion -Wno-unused-dummy-argument -fbacktrace -fcheck=all $(OMP)

#FOLDERS
OUT = ./output
BIN = ./bin
SRC = ./sources

#SOURCES
MAIN = $(SRC)/RTS_main.f90
FSOURCES = $(BIN)/RTS_global.o     $(BIN)/RTS_input.o $(BIN)/RTS_scattering.o $(BIN)/RTS_absorption.o \
	   $(BIN)/RTS_functions.o  $(BIN)/RTS_bc.o    $(BIN)/RTS_start.o      $(BIN)/RTS_solvers.o    \
	   $(BIN)/RTS_output.o     $(BIN)/RTS_radiation.o $(BIN)/RTS_energy.o 

all: presentation create $(PROGRAM)

$(PROGRAM): $(MAIN) $(FSOURCES)
	$(FC) $(OMP) $(FSOURCES) $(MAIN) -I $(BIN) -o $(BIN)/$(PROGRAM)
	@echo ' '
	@echo '--- Finished Building Program: $@ ---'

$(BIN)/%.o: $(SRC)/%.f90
	$(FC) $(FFLAGS) $< -J $(BIN) -o $@

run: create $(PROGRAM)
	./bin/$(PROGRAM)
	
clean:
	rm -f $(BIN)/*
	rm -f $(OUT)/*
	
delete:
	rm -r $(BIN)
	rm -r $(OUT)
	
create:
	@mkdir -p $(OUT) $(BIN)

cleanout:
	rm -f $(OUT)/*
	
#Backup
date=` data +%F `
version=` date +%d.%m.%y_%H.%M`

output      = ../$(PROGRAM)_$(version).tar.gz
saveindir   = ../$(PROGRAM)_INPUT_$(version).tar.gz
casedir     = ../$(PROGRAM)_CASE_$(version).tar.gz
files       = ./sources/*.f90 ./input/*.rts ./makefile
files_input = ./sources/RTS_functions.f90 ./input/*.rts
files_case  = ./sources/RTS_functions.f90 ./input/*.rts ./output/*
backup:
	@tar -czvf $(output) $(files)
	@echo '--- All program files have been saved ---'
	
saveinput: 
	@tar -czvf $(saveindir) $(files_input)
	@echo '--- All input files have been saved ---'
	
savecase: 
	@tar -czvf $(casedir) $(files_case)
	@echo '--- All case files have been saved ---'
	
help: presentation
	@echo '_________________________________________________________________________'
	@echo '                                                                         '
	@echo '                           Make Assistant Tool                           '
	@echo '_________________________________________________________________________'
	@echo '                               Make Flags                                '
	@echo '         make  < flag >                                                  '
	@echo '               run           : build and execute'
	@echo '               clean         : clears compilation files'
	@echo '               delete        : delete directories created'
	@echo '               cleanout      : clean results directory'
	@echo '               backup        : performs a code backup'
	@echo '               saveinput     : backup all input files'
	@echo '               savecase      : backup all case files'
	@echo '_________________________________________________________________________'

presentation:
ifeq ($(OS),Windows_NT)
	@echo '                                                                         '
	@echo '          ...........      ..................      ............          '
	@echo '          .............    ..................    ...............         '
	@echo '          ...       ....           ...          ....         ...         '
	@echo '          ...        ....          ...          ...                      '
	@echo '          ...        ....          ...          ....                     '
	@echo '          ...       ....           ...          ........                 '
	@echo '          ...    ......            ...           ............            '
	@echo '          ..........               ...                .........          '
	@echo '          ...    ....              ...                      .....        '
	@echo '          ...      ....            ...                       ....        '
	@echo '          ...       ....           ...         .             ....        '
	@echo '          ...        .....         ...         .....      ......         '
	@echo '          ...          ....        ...         ...............           '
	@echo '          ...           ....       ...             .........             '
else
	@echo '                                                                         '
	@echo '                   ███████████   ███████████  █████████                  '
	@echo '                  ░░███░░░░░███ ░█░░░███░░░█ ███░░░░░███                 '
	@echo '                   ░███    ░███ ░   ░███  ░ ░███    ░░░                  '
	@echo '                   ░██████████      ░███    ░░█████████                  '
	@echo '                   ░███░░░░░███     ░███     ░░░░░░░░███                 '
	@echo '                   ░███    ░███     ░███     ███    ░███                 '
	@echo '                   █████   █████    █████   ░░█████████                  '
	@echo '                  ░░░░░   ░░░░░    ░░░░░     ░░░░░░░░░                   '
endif
	@echo '                Radiative Transfer Simulator version 1.7.0               '
	@echo '                        Created by G. S. Rodrigues                       '
	@echo '                                                                         '

.PHONY: clean backup delete run
