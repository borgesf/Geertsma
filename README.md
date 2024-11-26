# Geertsma
Matlab implementation of Geertsma's equation (Geerstma, 1973) for subsidence and inflation of a reservoir in a homgeneous medium. In-line syntax is described in  the comments. See the files for further documentation.

- 'Geertsma_Checklist': runs a quick check on the functions.

- 'Geertsma_Exact': calculates the vertical compaction on the line above and below a compacting disc

- 'Geertsma_Example': runs a numerical example, with some plot suggestions. 

- 'Geertsma_Symbolic_ToolBox': modeling function, with detailed explanation. Use this if you have the symbolic toolbox installed

- 'Geertsma_No_ToolBox': modeling function, with detailed explanation. Use this if you do not have the symbolic toolbox installed. Accuracy is slightly worse.

- Other functions are numerical implementations necessary in the absence of the Symbolic Toolbox.

Author: 

	Filipe Borges 
    
Update History:

	- 14/11/2018: Added code to calculate Heuman Lambda Function (end of main code).
    - 23/11/2024: Update plot display

References:

    Fjær, E., R. M. Holt, A. Raaen, R. Risnes, and P. Horsrud,
       2008, Petroleum related rock mechanics: Elsevier, 53.

    Geertsma, J., 1973, A basic theory of subsidence due to
      reservoir compaction: the homogeneous case: Verhandelingen
      Kon. Ned. Geol. Mijnbouwk. Gen, 28, 43-62.

    https://medium.com/@FilipeBorgesBR/modelling-reservoir-induced-displacement-using-geertsmas-equation-in-matlab-db6decb3debc
