Set up:
  30-4-2025

Build with 
  Django
  MySQL
  Jmol/Jsmol

Database connection: 
  variantPSite/setting.py -> databases
  Mysql database in: msc.bc.ic.ac.uk
        Database name: variantp   User name: variantp   
        mysql -h msc.bc.ic.ac.uk  -u variantp -p variantp

Backend: 
  dbQuery/views.py 
  dbQuery/utils.py
  dbQuery/urls.py
  dbQuery/Class.py

Static: 
  dbQuery/static/dbQuery
  *Some AlphaFold models are uploaded here, but not all. 

Frontend: 
  \dbQuery\templates
  - main result page: /dbQuery/templates/dbQuery/results_2.html
  - Jsmol page: /dbQuery/templates/dbQuery/jsmol

