git checkout gh-pages
cd ..
javadoc -d gh-pages/ -sourcepath src net.sf.mendelscan
cd gh-pages
git add --all
git commit -m 'updated javadoc'
git push
cd ..
git commit -m 'update gh-pages submodule reference' gh-pages
git push
