get_diff(){
  # Execute the command and capture its output
  output=$("${CI_PROJECT_DIR}/bin/$1" $2)

  # Extract the Difference value
  value=$(echo "$output" | grep "Difference after" | cut -d: -f2)

  # Return the extracted value
  echo "$value"
}

compare(){
  if [ "$1" != "$2" ]; then
    echo "Comparing failed ($1 != $2)"
    exit 1
  fi
}

echo "Checking exemple execution..."

echo "RB Gauss-Seidel no PATSMA"
base=$(get_diff noat 1)

echo "RB Gauss-Seidel outside PATSMA"
echo "Dimension 1"
value=$(get_diff outsideat 1)
compare $value $base
echo "Dimension 2"
value=$(get_diff outsideat 2)
compare $value $base

echo "RB Gauss-Seidel inside PATSMA"
echo "Dimension 1"
value=$(get_diff insideat 1)
compare $value $base
echo "Dimension 2"
value=$(get_diff insideat 2)
compare $value $base

