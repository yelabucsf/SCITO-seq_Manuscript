'''
Clara clustering implementation
C code taken from https://github.com/cran/cluster
    Clustering LARge Applications
     ~		~~~   ~
     Clustering program based upon the k-medoid approach,
     and suitable for data sets of at least 100 objects.
     (for smaller data sets, please use program pam.)
'''
import clara
status = clara.system("ls -l")