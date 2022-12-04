read -p 'Enter filename: >> ' file
read -p 'Options: (1) Plot (2) Edit (3) Clean >> ' option
if [ $option -eq 2 ]
then
  echo -e 'Renaming file and table object ...'
  acol -read $file markov_chain_0 -set obj_name mc_0 -internal main.o2
  echo -e 'Success: Saved file as main.o2'
  read -p 'Actions: (1) Filter (2) Function (3) Print (4) Convert (5) Concatenate >> ' action
  if [ $action -eq 1 ]
  then
    read -p 'Filter by: (1) Rows (2) Columns >> ' filter
    if [ $filter -eq 1 ]
    then
      read -p 'Filter by Rows: (1) Zeros (2) N-rows (3) N-walker (4) Condition >> ' filter
      if [ $filter -eq 1 ]
      then
        echo -e 'Removing zero rows ...'
        acol -read tmp mc0 -select-rows "mult>0" -internal tmpnz
        echo -e 'Done.'
        read -p 'Filter by Rows: (1) Continue (2) Return (3) Exit >> ' choice
        if [ $choice -eq 1 ]
        then
          read -p 'Filter Rows: (1) N-rows (2) N-walker (3) Condition >> ' filter
          if [ $filter -eq 1 ]
          then
            read -p 'Filter by N-rows: Nmax = ' nindex
            echo -e 'Removing rows N<Nmax ...'
            acol -read tmpnz mc0 -index -delete-rows "N<$nindex" -internal tmpnr 
            echo -e 'Done.'
            read -p 'Filter Rows: (1) Continue (2) Return (3) Exit >> ' choice
            if [ $choice -eq 1 ]
            then
              read -p 'Filter Rows: (1) N-walker (2) Condition >> ' filter
              if [ $filter -eq 1 ]
              then
                read -p 'Filter Rows by N-walker: walker = ' nwalk
                echo -e 'Selecting rows walker==N ...'
                acol -read tmpnr mc0 -select-rows "walker==$nwalk" -internal tmpnw 
                echo -e 'Done.'
                read -p 'Filter by Rows: (1) Continue (2) Return (3) Exit >> ' choice
                if [ $choice -eq 1 ]
                then
                  read -p 'Filter Rows by condition = ' condition
                  echo -e 'Selecting rows by set condition ...'
                  acol -read tmpnw mc0 -select-rows "$condition" -internal tmpcr 
                  echo -e 'Done.'
                  read -p 'Filter by Rows: (1) Return (2) Exit >> ' choice
