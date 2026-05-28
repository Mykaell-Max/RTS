X_CO2 = 0.85
X_H2O = 0.10
X_N2 = 0.05


# X_CO2 = 1.0e-4
# X_H2O = 0.5
# X_N2 = 1 - X_CO2 - X_H2O

MM_CO2 = 44.009
MM_H2O = 18.015
MM_N2 = 28.014

denominator = MM_N2 * X_N2 + MM_CO2 * X_CO2 + MM_H2O * X_H2O

Y_CO2 = (X_CO2 * MM_CO2) / denominator
Y_H2O = (X_H2O * MM_H2O) / denominator
Y_N2 = (X_N2 * MM_N2) / denominator

print(Y_CO2, Y_H2O, Y_N2)
print(f"CO2:{Y_CO2}, H2O: {Y_H2O}, N2: {Y_N2}")