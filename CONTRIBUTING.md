Contributing to summit branch:

  1. click the "Fork" button on https://github.com/TorkamaniLab/Imputation_Autoencoder
  2. clone your local copy:

     cd /ccs/proj/bif138/$USER
     git clone -b summit git@github.com:<userid>/Imputation_Autoencoder
     git remote add upstream git@github.com:TorkamaniLab/Imputation_Autoencoder

  3. create a feature branch:

     git checkout -b <new_feature_branch_name>

  4. edit, test, then commit

     git add <changed file names>
     git commit -m 'description'

  5. check for changes to upstream

     git fetch
     git rebase summit # fix your commits to match any changes to summit

  5. push changes to your github

     git push -u origin <new_feature_branch_name>

  6. Use the pull request button on: https://github.com/<userid>/Imputation_Autoencoder

  7. change back to summit branch and continue

     git checkout summit
