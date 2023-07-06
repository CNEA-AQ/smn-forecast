#!/bin/bash

title="SILAMv5 User Manual" 
subtitle=`LC_ALL=en_US.UTF-8 date +"%a %d %b, %Y "` #date` # "10/10/2022"
version="5.0"

ouf="SILAMv${version}-UserGuide.pdf"
inf="silam_ug_v5.md"

pandoc -s -N --template=./templates/mytemplate.tex --variable mainfont="Liberation Sans" --variable sansfont="Liberation Sans Narrow" --variable monofont="Liberation Mono" --variable fontsize=14pt --variable version="${version}" -fmarkdown-implicit_figures --variable title="${title}" --variable subtitle="${subtitle}" --toc --variable geometry:margin=0.5in --pdf-engine=xelatex -s -o ${ouf} ${inf}

#pandoc -s -N --template=./templates/mytemplate.tex --filter ./filter/comments.py --variable mainfont="Times New Roman" --variable sansfont="Helvetica" --variable monofont="Menlo" --variable fontsize=12pt --variable version=5.4 -fmarkdown-implicit_figures --variable title="CMAQv5.4 User Manual" --variable subtitle="10/10/2022" --toc --variable geometry:margin=1in --pdf-engine=xelatex -s -o ./PDF/SILAM_UG.pdf silam_ug_v5.md

