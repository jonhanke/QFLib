# Introduction #

Below are instructions for using the Google code mercurial repository.  For general Mercurial instructions, you can look at the excellent online book

> http://hgbook.red-bean.com/read/mercurial-in-daily-use.html

or the quickstart guide

> http://mercurial.selenic.com/quickstart/

on the official mercurial website.


# Instructions for the Google Code  Repository #

To clone the Google Code repository, use:

> hg clone https://code.google.com/p/qflib/

To push changes to the Google Code repository, enter the repository directory and use:

> hg push https://code.google.com/p/qflib/





# Instructions for the Sourceforge Repository (not maintained) #

Due to some initial trouble pushing changes to Google code using Mercurial, an alternate mercurial repository was created at Sourceforge.net.  Now that the Google code repository is working, this one is not maintained.


To clone the SourceForge repository, use:

> hg clone ssh://hg.code.sf.net/p/qflibrary/code qflibrary-code


To push changes to the SourceForge repository, enter the repository directory and use:

> hg push ssh://hg.code.sf.net/p/qflibrary/code