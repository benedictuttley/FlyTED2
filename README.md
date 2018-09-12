# FlyTED2 - Drosophila Testes Expression Data Base

A database holding Drosophila (fruit fly) testes gene expression images and annotations.
A project to allow simple retrieval of image data and the images respective annotations for > 1000 genes, the annotations are
supported by additional data pulled in from FlyMine using the Intermine API: http://intermine.readthedocs.io/en/latest/web-services

Users can search for genes and provide an optional variant argument (which mutant of that gene you are interested in), the application then returns:

- Transcripts and probe sequences
- Images and metadata
- Descriptions
- Graph showing expression by tissue: Credit to Flymine for this data set.
- Graph showing expression by variant

All code is released under the MIT license, the actual expression images themselves are not held on GitHub due to the large quantity, however, they are published under the Attribution 2.0 UK: England & Wales, so if to be used elsewhere, then the individual needs to credit the authors listed on the site.
