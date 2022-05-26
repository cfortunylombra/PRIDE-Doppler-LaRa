#!/bin/bash - screen+bash - run code with source
echo "Activate tudat-bundle env"
conda activate tudat-bundle
# Reference RISE+LaRa+PRIDE+corr=0

# # Case 1 (only RISE)
# echo "---CASE 1 (only RISE)---"
# echo "Copy py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE.py"
# cp -r ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE.py ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo "Change build path"
# sed -i '14s+sys.path.insert(0, "/home/cfortunylombra/tudat-bundle/cmake-build-release-wsl/tudatpy/")+sys.path.insert(0, "/home/carlos/tudat-bundle/build_release/tudatpy/")+' ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo "LaRa_boolean = True -> False"
# sed -i '42s+LaRa_boolean = True+LaRa_boolean = False+' ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo "Run code"
# python ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo ""

# # Case 2 (only LaRa)
# echo "---CASE 2 (only LaRa PRIDE-False remove-False - corr=0)---"
# echo "Copy py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE.py"
# cp -r ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE.py ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo "Change build path"
# sed -i '14s+sys.path.insert(0, "/home/cfortunylombra/tudat-bundle/cmake-build-release-wsl/tudatpy/")+sys.path.insert(0, "/home/carlos/tudat-bundle/build_release/tudatpy/")+' ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo "RISE_boolean = False -> True"
# sed -i '41s+RISE_boolean = True+RISE_boolean = False+' ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo "PRIDE_boolean = True -> False"
# sed -i '46s+PRIDE_boolean = True+PRIDE_boolean = False+' ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo "Run code"
# python ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo ""

# # Case 3 (only LaRa)
# echo "---CASE 3 (only LaRa PRIDE-True remove-True - corr=0)---"
# echo "Copy py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE.py"
# cp -r ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE.py ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo "Change build path"
# sed -i '14s+sys.path.insert(0, "/home/cfortunylombra/tudat-bundle/cmake-build-release-wsl/tudatpy/")+sys.path.insert(0, "/home/carlos/tudat-bundle/build_release/tudatpy/")+' ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo "RISE_boolean = False -> True"
# sed -i '41s+RISE_boolean = True+RISE_boolean = False+' ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo "remove_PRIDE_weight_boolean = False -> True"
# sed -i '49s+remove_PRIDE_weight_boolean = False+remove_PRIDE_weight_boolean = True+' ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo "Run code"
# python ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo ""

# # Case 4 (only LaRa)
# echo "---CASE 4 (only LaRa PRIDE-True remove-False - corr=0)---"
# echo "Copy py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE.py"
# cp -r ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE.py ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo "Change build path"
# sed -i '14s+sys.path.insert(0, "/home/cfortunylombra/tudat-bundle/cmake-build-release-wsl/tudatpy/")+sys.path.insert(0, "/home/carlos/tudat-bundle/build_release/tudatpy/")+' ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo "RISE_boolean = False -> True"
# sed -i '41s+RISE_boolean = True+RISE_boolean = False+' ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo "Run code"
# python ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo ""

# # Case 5 (only LaRa)
# echo "---CASE 5 (only LaRa PRIDE-True remove-False - corr=0.99)---"
# echo "Copy py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE.py"
# cp -r ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE.py ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo "Change build path"
# sed -i '14s+sys.path.insert(0, "/home/cfortunylombra/tudat-bundle/cmake-build-release-wsl/tudatpy/")+sys.path.insert(0, "/home/carlos/tudat-bundle/build_release/tudatpy/")+' ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo "RISE_boolean = False -> True"
# sed -i '41s+RISE_boolean = True+RISE_boolean = False+' ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo "correlation = 0 -> 0.99"
# sed -i '52s+correlation = 0+correlation = 0.99+' ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo "Run code"
# python ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo ""

# # Case 7 (RISE + LaRa)
# echo "---CASE 7 (RISE+LaRa PRIDE-False remove-False - corr=0)---"
# echo "Copy py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE.py"
# cp -r ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE.py ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo "Change build path"
# sed -i '14s+sys.path.insert(0, "/home/cfortunylombra/tudat-bundle/cmake-build-release-wsl/tudatpy/")+sys.path.insert(0, "/home/carlos/tudat-bundle/build_release/tudatpy/")+' ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo "PRIDE_boolean = True -> False"
# sed -i '46s+PRIDE_boolean = True+PRIDE_boolean = False+' ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo "Run code"
# python ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo ""

# # Case 8 (RISE + LaRa)
# echo "---CASE 8 (RISE+LaRa PRIDE-True remove-True - corr=0)---"
# echo "Copy py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE.py"
# cp -r ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE.py ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo "Change build path"
# sed -i '14s+sys.path.insert(0, "/home/cfortunylombra/tudat-bundle/cmake-build-release-wsl/tudatpy/")+sys.path.insert(0, "/home/carlos/tudat-bundle/build_release/tudatpy/")+' ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo "remove_PRIDE_weight_boolean = False -> True"
# sed -i '49s+remove_PRIDE_weight_boolean = False+remove_PRIDE_weight_boolean = True+' ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo "Run code"
# python ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo ""

# # Case 9 (RISE + LaRa)
# echo "---CASE 9 (RISE+LaRa PRIDE-True remove-False - corr=0)---"
# echo "Copy py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE.py"
# cp -r ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE.py ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo "Change build path"
# sed -i '14s+sys.path.insert(0, "/home/cfortunylombra/tudat-bundle/cmake-build-release-wsl/tudatpy/")+sys.path.insert(0, "/home/carlos/tudat-bundle/build_release/tudatpy/")+' ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo "Run code"
# python ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo ""

# # Case 10 (RISE + LaRa)
# echo "---CASE 10 (RISE+LaRa PRIDE-True remove-False - corr=0.99)---"
# echo "Copy py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE.py"
# cp -r ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE.py ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo "Change build path"
# sed -i '14s+sys.path.insert(0, "/home/cfortunylombra/tudat-bundle/cmake-build-release-wsl/tudatpy/")+sys.path.insert(0, "/home/carlos/tudat-bundle/build_release/tudatpy/")+' ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo "correlation = 0 -> 0.99"
# sed -i '52s+correlation = 0+correlation = 0.99+' ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo "Run code"
# python ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo ""

for var_corr in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
do 
	# Case 6 (only LaRa)
	echo "---CASE 6 (only LaRa PRIDE-True remove-False - corr=$var_corr)---"
	echo "Copy py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE.py"
	cp -r ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE.py ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
	echo "Change build path"
	sed -i '14s+sys.path.insert(0, "/home/cfortunylombra/tudat-bundle/cmake-build-release-wsl/tudatpy/")+sys.path.insert(0, "/home/carlos/tudat-bundle/build_release/tudatpy/")+' ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
	echo "RISE_boolean = False -> True"
	sed -i '41s+RISE_boolean = True+RISE_boolean = False+' ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
	echo "correlation = 0 -> $var_corr"
	sed -i "52s+correlation = 0+correlation = $var_corr+" ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
	echo "Run code"
	python ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
	echo ""
done

for var_corr in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
do 
	# Case 11 (RISE + LaRa)
	echo "---CASE 11 (RISE+LaRa PRIDE-True remove-False - corr=$var_corr)---"
	echo "Copy py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE.py"
	cp -r ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE.py ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
	echo "Change build path"
	sed -i '14s+sys.path.insert(0, "/home/cfortunylombra/tudat-bundle/cmake-build-release-wsl/tudatpy/")+sys.path.insert(0, "/home/carlos/tudat-bundle/build_release/tudatpy/")+' ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
	echo "correlation = 0 -> $var_corr"
	sed -i "52s+correlation = 0+correlation = $var_corr+" ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
	echo "Run code"
	lspython ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
	echo ""
done