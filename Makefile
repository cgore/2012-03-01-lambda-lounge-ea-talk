NAME=slides

all: $(NAME).dvi $(NAME).ps $(NAME).pdf

$(NAME).dvi: $(NAME).latex
	latex $(NAME).latex

$(NAME).ps: $(NAME).dvi
	dvips -o $(NAME).ps -t landscape $(NAME).dvi

$(NAME).pdf: $(NAME).ps
	ps2pdf $(NAME).ps

print: $(NAME).ps
	lpr $(NAME).ps
