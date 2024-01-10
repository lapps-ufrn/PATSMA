echo "Checking code integrity..."
${CI_PROJECT_DIR}/bin/integrity

echo "Checking exemple execution..."
echo "RB Gauss-Seidel outside PATSMA"
${CI_PROJECT_DIR}/bin/outsideat 1
${CI_PROJECT_DIR}/bin/outsideat 2
echo "RB Gauss-Seidel inside PATSMA"
${CI_PROJECT_DIR}/bin/insideat 1
${CI_PROJECT_DIR}/bin/insideat 2
echo "RB Gauss-Seidel no PATSMA"
${CI_PROJECT_DIR}/bin/noat