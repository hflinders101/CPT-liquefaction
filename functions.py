import pandas as pd
import numpy as np
import warnings
import scipy.integrate as integrate

def soil_parameters(df):
    Pa = 101.325  # Atmospheric pressure in kPa

    # /////////////////////////////////////////////// COLUMNS \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Add new columns for each of the soil parameters
    new_columns = ['u calc','qc calc', 'qt calc', 'Qt', "Rf (%)", "Gamma (kN/m^3)", "Total Stress (kPa)",
                   "Effective Stress (kPa)",
                   "Fr (%)", 'n1', "Cn", "Qtn", "Ic", 'n2', 'error', 'OCR R', 'OCR K', "cu_bq", "cu_14", "M", "k0_1",
                   'k0_2', "Vs R",
                   'Vs M', "k (m/s)", 'ψ', "φ' R", "φ' K", "φ' J", 'Qtn,cs', "φ' M", "φ' U", 'Dr B', 'Dr K', 'Dr J',
                   'Dr I',
                   'Cn2', "qc1", 'qc2', 'error2']
    df_new_columns = pd.DataFrame(columns=new_columns)
    df = pd.concat([df, df_new_columns], axis=1)

    # /////////////////////////////////////////////// end COLUMNS \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

    # ///////////////////////////////////////////// GENERAL CALCULATIONS \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # "qc calc" and "qt calc" columns created to change bad data into numbers our equations can handle. Units are also converted
    df['qc calc'] = df['qc (MPa)'] * 1000
    df['qt calc'] = df['qt (MPa)'] * 1000
    for i in range(len(df.index)):
        row = df.loc[i].copy(deep=False)
        if row['qc calc'] <= 0:
            df.at[i, 'qc calc'] = float('NaN')
        if row['qt calc'] <= 0:
            df.at[i, 'qt calc'] = float('NaN')

    # Rf calc
    def calcRf(fs, qt_calc):
        if fs < 0.00001:
            return 0
        elif np.isnan(qt_calc):
            return 0
        else:
            return np.divide(fs, qt_calc) * 100

    df['Rf (%)'] = [calcRf(x, y) for x, y in zip(df['fs (kPa)'], df['qt calc'])]

    # Gamma calc
    def calcGamma(Rf, qt_calc):
        if Rf <= 0:
            return 18.08  # default gamma value when there's a pre-hole
        else:
            return 9.81 * (0.27 * np.log10(Rf) + 0.36 * np.log10(qt_calc / Pa) + 1.236)

    df['Gamma (kN/m^3)'] = [calcGamma(x, y) for x, y in zip(df['Rf (%)'], df["qt calc"])]

    # Total Stress calculation
    df['Total Stress (kPa)'] = df['Gamma (kN/m^3)'] * df['Depth (m)']
    for i in range(1, len(df.index)):
        row = df.loc[i]
        previous = df.loc[i - 1]
        df.at[i, 'Total Stress (kPa)'] = (row['Depth (m)'] - previous['Depth (m)']) * row['Gamma (kN/m^3)'] + previous[
            'Total Stress (kPa)']

    # Effective Stress calculation
    if df.loc[0]['GWT [m]'] > 0:
        GWT = df.loc[0]['GWT [m]']
        df['Effective Stress (kPa)'] = df['Total Stress (kPa)']

        df['u calc'] = 0

        for i in range(len(df.index)):
            row = df.loc[i]
            if row['Depth (m)'] >= GWT:
                df.at[i, 'Effective Stress (kPa)'] = row['Total Stress (kPa)'] - ((row['Depth (m)'] - GWT) * 9.81)
                df.at[i, 'u calc'] = ((row['Depth (m)'] - GWT) * 9.81).astype(np.int64)
                # print(type(row["u (kPa)"]), type(row["u (kPa)"]))
            # Fr calcuation
            if row['fs (kPa)'] <= 0:
                df.at[i, 'Fr (%)'] = 0
            else:
                df.at[i, 'Fr (%)'] = (row["fs (kPa)"] / (row["qt calc"] - row['Total Stress (kPa)']) * 100).astype(
                    float)
    else:
        warnings.warn('GWT marked as 0 or not provided')

    # Qt calculation
    df['Qt'] = [(x - y) / z for x, y, z in zip(df["qt calc"], df["Total Stress (kPa)"], df['Effective Stress (kPa)'])]
    # ///////////////////////////////////////////// end GENERAL CALCULATIONS \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

    # ////////////////////////////////////////////// Ic CALCULATION \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    df['n1'] = 1  # Use 1 as the first guess for n
    tolerance = 0.01  # Define the Ic iteration tolerance here
    counter = False

    while not counter:
        # Calculate Cn
        df['Cn'] = (Pa / df['Effective Stress (kPa)']) ** df['n1']
        df['Cn'] = [1.7 if x >= 1.7 else x for x in df['Cn']]

        # Calculate Qtn
        df['Qtn'] = (((df["qt calc"] - df['Total Stress (kPa)']) / Pa) * df['Cn']).astype(float)

        # Calculate Ic
        for i in range(len(df.index)):
            row = df.loc[i]
            if row['Fr (%)'] <= 0 or row['Qtn'] <= 0:
                df.at[i, 'Ic'] = 0
            elif row['Fr (%)'] == float('NaN') or row['Qtn'] == float('NaN'):
                df.at[i, 'Ic'] = float('NaN')
            else:
                df.at[i, 'Ic'] = (((3.47 - np.log10(row['Qtn'])) ** 2) + (np.log10(row['Fr (%)']) + 1.22) ** 2) ** 0.5

        # Calculate n2
        for i in range(len(df.index)):
            row = df.loc[i]
            temp = 0.381 * (row['Ic']) + 0.05 * (row['Effective Stress (kPa)'] / Pa) - .15
            if temp > 1:
                df.at[i, 'n2'] = 1
            else:
                df.at[i, 'n2'] = temp

        # Calculate the error and set n2 as n1 for further iterations
        df['error'] = df['n1'] - df['n2']
        df['n1'] = df['n2']

        # Check to see if every row meets our error tolerance. If not, repeat the process
        counter1 = True
        for i in range(len(df.index)):
            row = df.loc[i]
            if row['Ic'] > 0 and row['error'] > tolerance:
                counter1 = False
                break
        counter = counter1
    # /////////////////////////////////////////// end Ic CALCULATION \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

    # //////////////////////////////// Dr CALCULATION Idriss and Boulanger 2008 \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

    df['qc1'] = df['qc calc']  # Set recorded qc values as initial qc1n guess
    tolerance = 0.01  # Define the Dr iteration tolerance here

    counter = False
    it_counter = 0  # Create variable to count the number of iterations

    while not counter:
        # Cn calculation
        df['Cn2'] = (Pa / df['Effective Stress (kPa)']) ** (1.338 - .249 * df['qc1'] ** .264)

        # New qcn1 calculation
        df['qc2'] = df['Cn2'] * df['qc calc'] / Pa

        # Dr calculation
        df['Dr I'] = .478 * df['qc1'] ** .264 - 1.063

        # Find error between guess and new qcn1 calculation
        df['error2'] = np.abs(df['qc1'] - df['qc2'])

        # Set new qc1n calculation as new guess for next iteration
        df['qc1'] = df['qc2']

        # Count the number of iterations
        it_counter += 1

        # Check to see if every row meets our error tolerance. If not, repeat the process.
        # If there have been more than 100 iterations, set the value to "No Solution"
        counter1 = True
        for i in range(len(df.index)):
            row = df.loc[i]
            if it_counter == 100:
                if row['Dr I'] > 0 and row['error2'] > tolerance:
                    df.at[i, 'Dr I'] = str('No Solution')
            else:
                if row['Dr I'] > 0 and row['error2'] > tolerance:
                    counter1 = False
                    break
        counter = counter1

    # //////////////////////////////////// end Dr CALCULATION Idriss and Boulanger 2008 \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

    # //////////////////////////////////////////// COHESIVE LAYER PROPERTIES \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

    # ---------------------------------------- OCR calculations --------------------------------------------------------
    # Robertson 2009
    def calcOCR_R(Ic, Qt):
        if Ic >= 2.6:
            return .25 * Qt ** 1.25

    df['OCR R'] = [calcOCR_R(x, y) for x, y in zip(df['Ic'], df['Qt'])]

    # Kulkawy and Mayne 1990
    def calcOCR_K(Ic, Qt):
        k = 0.33  # An average value of k = 0.33 can be assumed, with an expected range of 0.2 to 0.5. Higher values of k are recommended in aged, heavily overconsolidated clays.
        if Ic >= 2.6 and Qt < 20:
            return k * Qt

    df['OCR K'] = [calcOCR_K(x, y) for x, y in zip(df['Ic'], df['Qt'])]
    # -----------------------------------end OCR calculations ----------------------------------------------------------

    # Begin for loop to perform cell based calculations
    for i in range(len(df.index)):
        row = df.loc[i]
        if row['Ic'] >= 2.6:
            df.at[i, 'Dr I'] = float('NaN')
        if row['Ic'] == 0:
            df.at[i, 'Dr I'] = float('NaN')

        if row["Ic"] >= 2.6:  # Check the soil type

            # --------------------------- cu calculations --------------------------------------------------------------
            # Mayne & Peuchen 2018
            GWT = df.loc[0]['GWT [m]']
            if row['Depth (m)'] >= GWT:
                u0 = (row['Depth (m)'] - GWT) * 9.81
            else:
                u0 = 0
            Bq = (row['u calc'] - u0) / (row['qt calc'] - row['Total Stress (kPa)'])
            if Bq <= -0.1:
                Bq = -0.009999999
            Nkt = 10.5 - 4.6 * np.log(Bq + 0.1)
            df.at[i, 'cu_bq'] = (row['qt calc'] - row['Total Stress (kPa)']) / Nkt
            df.at[i, 'cu_14'] = (row['qt calc'] - row[
                'Total Stress (kPa)']) / 14  # Dr. Rollins wanted to use a set value of Nkt = 14 in addition to the bq calc since he is unfamiliar with bq
            # -------------------------- end cu calculations -----------------------------------------------------------

            # ----------------------------- M calculations -------------------------------------------------------------
            # Robertson 2009. From what I can tell from the paper, M is in MPa
            if row['Qt'] >= 14:
                df.at[i, 'M'] = (row['qt calc'] - row['Total Stress (kPa)']) * 14
            else:
                df.at[i, 'M'] = (row['qt calc'] - row['Total Stress (kPa)']) * row['Qt']
            # ------------------------------- end M calculations -------------------------------------------------------

            # -------------------------------k0 calculations -----------------------------------------------------------
            # Kulhway and Mayne 1990
            df.at[i, 'k0_1'] = (row['qt calc'] - row['Total Stress (kPa)']) / row['Effective Stress (kPa)'] * .1
            df.at[i, 'k0_2'] = 0.5 * (row['OCR R']) ** 0.5
            # -------------------------------end k0 calculations -------------------------------------------------------

            # ------------------------------- Vs calculation -----------------------------------------------------------
            # Robertson 2009
            avs = 10 ** (0.55 * row['Ic'] + 1.68)
            if (avs * (row['qt calc'] - row['Total Stress (kPa)'])) > 0:
                df.at[i, 'Vs R'] = (avs * (row['qt calc'] - row['Total Stress (kPa)']) / Pa) ** .5

            # Mayne 2006
            if row['fs (kPa)'] > 0:
                df.at[i, 'Vs M'] = 51.6 * np.log(row['fs (kPa)']) + 18.5
            # ---------------------------------end Vs calculation ------------------------------------------------------

            # --------------------------------k for permeability -------------------------------------------------------
            # Robertson 2015
            if row['Ic'] < 3.27:
                df.at[i, 'k (m/s)'] = 10 ** (.952 - 3.04 * row['Ic'])
            if 3.27 < row['Ic'] < 4:
                df.at[i, 'k (m/s)'] = 10 ** (-4.52 - 1.37 * row['Ic'])
            # --------------------------------end k for permeability ---------------------------------------------------

            # ------------------------------- φ' calculation -----------------------------------------------------------
            # Mayne 2006
            GWT = df.loc[0]['GWT [m]']
            if row['Depth (m)'] >= GWT:
                u0 = (row['Depth (m)'] - GWT) * 9.81
            else:
                u0 = 0
            Bq = (row['u calc'] - u0) / (row['qt calc'] - row['Total Stress (kPa)'])
            if Bq <= 0:
                Bq = 0.1
            elif Bq > 1:
                Bq = 1
            if row['Qt'] > 0:
                df.at[i, "φ' M"] = 29.5 * Bq ** 0.121 * (0.256 + 0.336 * Bq + np.log10(row['Qt']))
            # ----------------------------- end φ' calculation ---------------------------------------------------------
    # /////////////////////////////////////// end COHESIVE LAYER PROPERTIES \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

    # /////////////////////////////////////// NON-COHESIVE LAYER PROPERTIES \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Begin for loop to perform cell based calculations
    for i in range(len(df.index)):
        row = df.loc[i]
        if 2.6 > row['Ic'] > 0:  # Check the soil type. Ic == 0 means there's not data.

            # ---------------------------------------- φ' calculation --------------------------------------------------
            # Robertson and Campanella 1983
            if row['qc calc'] > 0:
                df.at[i, "φ' R"] = np.degrees(
                    np.arctan(1 / 2.68 * (np.log10(row['qc calc'] / row['Effective Stress (kPa)']) + 0.29)))

            # Kulhawy and Mayne 1990
            df.at[i, "φ' K"] = 17.6 + 11 * np.log10(row['Qtn'])

            # Jefferies and Been 2006
            if row['Ic'] <= 1.64:
                Kc = 1.0
            elif 1.64 < row['Ic'] < 2.36 and row['Fr (%)'] < 0.5:
                Kc = 1.0
            elif 1.64 < row['Ic'] <= 2.5:
                Kc = 5.58 * (row["Ic"]) ** 3 - 0.403 * (row["Ic"]) ** 4 - 21.63 * (row["Ic"]) ** 2 + 33.75 * (
                    row["Ic"]) - 17.88
            else:
                Kc = 6 * 10 ** -7 * row['Ic'] ** 16.76
            df.at[i, "φ' J"] = 33 + 15.84 * (
                np.log10(Kc * row['Qtn'])) - 26.88  # Used a φ'cv value of 33 degrees per Dr. Rollins' instructions

            # df.at[i, 'Qtn,cs'] = Kc * row['Qtn'] -------- if we want to check Qtn,cs values for checking here you go

            # Uzielli, Mayne, and Cassidy 2013
            df.at[i, "φ' U"] = 25 * (row['qt calc'] / (row['Effective Stress (kPa)']) ** 0.5) ** 0.1
            # ------------------------------------- end φ' calculation -------------------------------------------------

            # ------------------------------------------- DR calculation -----------------------------------------------
            # Baldi et al. 1986      ******WEIRD NUMBERS**********
            C0, C2 = 15.7, 2.41  # For moderately compressible, normally consolidated, unaged and uncemented, predominantly quartz sands the constants are: C0 = 15.7 and C2 = 2.41
            Qcn = (row['qc calc'] / Pa) / (row['Effective Stress (kPa)'] / Pa) ** 0.5
            df.at[i, 'Dr B'] = (1 / C2) * np.log(Qcn / C0)

            # Kulhawy and Mayne 1990
            df.at[i, 'Dr K'] = (row['Qtn'] / 350) ** 0.5  # Used the Qtn/350 simplification of this equation per
            # Dr. Rollins' instructions since we don't have the needed
            # information for the non-simplified version of the equation

            # Jamiolkowski et al. 2003
            c0 = 17.68
            c1 = 0.5
            c2 = 3.10
            df.at[i, 'Dr J'] = 1 / c2 * np.log(
                (row['qt calc'] / Pa) / (c0 * (row['Effective Stress (kPa)'] / Pa) ** c1))
            # -------------------------------------- end DR calculation ------------------------------------------------

            # ---------------------------------- ψ state parameter calculation -----------------------------------------
            # Robertson 2010
            df.at[i, 'ψ'] = 0.56 - 0.33 * np.log10(Kc * row['Qtn'])
            # ---------------------------------- end  ψ state parameter calculation ------------------------------------

            # ------------------------------------- Vs calculation -----------------------------------------------------
            # Robertson 2009
            avs = 10 ** (0.55 * row['Ic'] + 1.68)
            df.at[i, 'Vs R'] = (avs * (row['qt calc'] - row['Total Stress (kPa)']) / Pa) ** 0.5

            # Mayne 2006
            if row['fs (kPa)'] > 0:
                df.at[i, 'Vs M'] = 51.6 * np.log(row['fs (kPa)']) + 18.5
            # ------------------------------------- end Vs calculation -------------------------------------------------

            # --------------------------------- k for permeability -----------------------------------------------------
            # Robertson 2010
            df.at[i, 'k (m/s)'] = 10 ** (0.952 - 3.04 * row['Ic'])
            # --------------------------------- end k for permeability -------------------------------------------------

            # ----------------------------------------- M --------------------------------------------------------------
            # Robertson 2009. From what I can tell from the paper, M is in MPa
            if row['Ic'] > 2.2:
                if row['Qt'] >= 14:
                    df.at[i, 'M'] = (row['qt calc'] - row['Total Stress (kPa)']) * 14
                else:
                    df.at[i, 'M'] = (row['qt calc'] - row['Total Stress (kPa)']) * row['Qt']
            else:
                am = 0.0188 * (10 ** (0.55 * row['Ic'] + 1.68))
                df.at[i, 'M'] = am * (row['qt calc'] - row['Total Stress (kPa)'])
            # ------------------------------------- end M --------------------------------------------------------------
        elif row['Ic'] == 0:
            df.at[i, 'Ic'] = float('NaN')

    # //////////////////////////////////////// end NON-COHESIVE LAYER PROPERTIES \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

    # Delete columns that stored variables for calculations but that we don't want in the final spreadsheet
    df.drop(['qc calc', 'qt calc', 'Qt', 'n1', 'Cn', 'Qtn', 'n2', 'error', 'qc1', 'qc2', 'error2', 'Cn2', 'Qtn,cs'],
            axis=1, inplace=True)
    return df

def PGA_insertion(df,PGA_filepath, site):

    pga = pd.read_excel(PGA_filepath)
    pga.set_index('site',inplace=True)
    df.at[0,'PGA_20may'] = pga.loc[site]['PGA_20may']
    df.at[0,'PGA_29may'] = pga.loc[site]['PGA_29may']
    df.at[0,'Liquefaction'] = pga.loc[site]['Liquefaction']
    return df

# input df must have PGA and Liquefaction values already defined
def FS_liq(df, Magnitude1, Magnitude2, date1, date2): # FS equation from Idriss and Boulanger 2008
    Pa = 101.325
    new_columns = ['qc1n','qc1ncs', 'Kσ', 'rd_'+date1, 'rd_'+date2, "CSR_"+date1, "CRR_"+date1, 'CSR_'+date2,
                   'CRR_'+date2, "FS_"+date1, "FS_"+date2]
    df_new_columns = pd.DataFrame(columns=new_columns)
    df = pd.concat([df, df_new_columns], axis=1)

    if df.loc[0]["GWT [m]"] < df.loc[0]['preforo [m]']:
        df.at[1, 'preforo [m]'] = 'preforo is below GWT'

    # FSliq part
    MSF1 = 6.9 * np.exp(-Magnitude1 / 4) - .058
    if MSF1 > 1.8:
        MSF1 = 1.8
    MSF2 = 6.9 * np.exp(-Magnitude2 / 4) - .058
    if MSF2 > 1.8:
        MSF2 = 1.8

    # Calculating K sigma
    for i in range(len(df.index)):
        row = df.loc[i]  # this takes a screenshot
        if row['Dr I'] == 'No Solution':
            df.at[i, 'qc1n'] = float('NaN')
        else:
            df.at[i, 'qc1n'] = ((row['Dr I'] + 1.063) / .478) ** (1 / .264)  # from Dr I iterative calc (we backcalculate here)
        row = df.loc[i]

        if 2.6 > row["Ic"] > 0:

            c_sigma = 1 / (37.3 - 8.27 * row['qc1n'] ** .264)
            if c_sigma > .3:
                c_sigma = .3

            Kσ = 1 - c_sigma * np.log(row["Effective Stress (kPa)"] / Pa)
            if Kσ > 1.1:
                Kσ = 1.1
            df.at[i, 'Kσ'] = Kσ

            # Calculating rd
            alpha = -1.012 - 1.126 * np.sin(
                row['Depth (m)'] / 11.73 + 5.133)  # rd is only good for depths less than 20 meters (pg 68)
            beta = .106 + .118 * np.sin(row['Depth (m)'] / 11.28 + 5.142)
            if row['Depth (m)'] < 20:
                df.at[i, 'rd_'+date1] = np.exp(alpha + beta * Magnitude1)
                df.at[i, 'rd_'+date2] = np.exp(alpha + beta * Magnitude2)

            row = df.loc[i]

            # Calcuating CSR
            g = 1
            df.at[i, "CSR_"+date1] = .65 * df.loc[0, "PGA_"+date1] / g * row["Total Stress (kPa)"] / row[
                "Effective Stress (kPa)"] * row["rd_"+date1] / MSF1 / row['Kσ']
            df.at[i, "CSR_"+date2] = .65 * df.loc[0, "PGA_"+date2] / g * row["Total Stress (kPa)"] / row[
                "Effective Stress (kPa)"] * row["rd_"+date2] / MSF2 / row['Kσ']

            row = df.loc[i]

            # Calcuatig CRR
            FC = 2 * 2.8 * row["Ic"]**2.6 #Taken from Emilia Romagna paper
            if FC > 100:
                FC = 100
            elif FC < 0:
                FC = 0
            qc1ncs = row["qc1n"] + (5.4 + row['qc1n'] / 16) * np.exp(1.63 + 9.7 / (FC + 0.01) - (15.7 / (FC + 0.01)) ** 2)
            df.at[i, "qc1ncs"] = qc1ncs
            # print(MSF1,row["Kσ"],i)
            df.at[i, "CRR_"+date1] = np.exp(
                qc1ncs / 540 + (qc1ncs / 67) ** 2 - (qc1ncs / 80) ** 3 + (qc1ncs / 114) ** 4 - 3) / MSF1 / row["Kσ"]

            df.at[i, "CRR_"+date2] = np.exp(
                qc1ncs / 540 + (qc1ncs / 67) ** 2 - (qc1ncs / 80) ** 3 + (qc1ncs / 114) ** 4 - 3) / MSF2 / row[
                                        "Kσ"]

            row = df.loc[i]

            # FS liq
            if row["Depth (m)"] <= df.loc[0]['GWT [m]']:
                df.at[i, "FS_"+date1] = 9999
                df.at[i, "FS_"+date2] = 9999
            else:
                df.at[i, "FS_"+date1] = row['CRR_'+date1] / row['CSR_'+date1]
                df.at[i, "FS_"+date2] = row['CRR_'+date2] / row['CSR_'+date2]

    return df


# calculates h2 as the thickness of the shallowest liquefiable layer greater than 0.3 meters
def h1_h2_basic (df, depth_column_name, FS_column_name):
  # Initialize variables
  last_liq_depth = None
  start_liq_depth = None
  h2_thickness = 0
  h1_thickness = df.iloc[-1][depth_column_name]


  # Iterate through the DataFrame
  for index, row in df.iterrows():
      depth = row[depth_column_name]
      FS = row[FS_column_name]
      if FS == '':
          FS = float('NaN')

      if FS < 1 and (last_liq_depth is None or depth - last_liq_depth <= 0.3):
          last_liq_depth = depth
          if start_liq_depth is None:
            start_liq_depth = depth
            if index == 0:
              h1_index = 0
            else:
              h1_index = index - 1
      else:
          if last_liq_depth is not None and depth - last_liq_depth > 0.3:
            h2_thickness = last_liq_depth - start_liq_depth
            h1_thickness = df.loc[h1_index][depth_column_name]
            break

  h1_column_name = "h1_basic" + FS_column_name.lstrip("FS")
  h2_columnn_name = "h2_basic" + FS_column_name.lstrip("FS")
  df.at[0,h1_column_name] = h1_thickness
  df.at[0,h2_columnn_name] = h2_thickness

  return df

# calculates h2 as the summation of all liquefiable layers for depths less than 10 meters
def h1_h2_cumulative(df, depth_column_name, FS_column_name):
    # Initialize variables
    current_depth = None
    start_depth = None
    h2_thickness = 0
    h1_thickness = 10

    # Iterate through the DataFrame
    for index, row in df.iterrows():
        depth = row[depth_column_name]
        FS = row[FS_column_name]
        if FS == '':
            FS = float('NaN')

        if FS < 1 and (current_depth is None or depth - current_depth <= 0.3):  # replace current_depth with start_depth
            # FS value found and consistent with the previous depth
            current_depth = depth
            if start_depth is None:
                start_depth = depth
                if index == 0:
                    h1_index = 0
                else:
                    h1_index = index - 1
        else:
            # FS value is not present or not consistent
            if current_depth is not None and current_depth - start_depth > 0.3:
                # Store the thickness as 'H1' in another part of the DataFrame
                h1_thickness = df.loc[h1_index][depth_column_name]
                break

            # Reset variables for the next consistent layer
            current_depth = None
            start_depth = None

    for index, row in df.iterrows():
        depth = row[depth_column_name]
        FS = row[FS_column_name]
        if FS == '':
            FS = float('NaN')
        if 0 < FS < 1 and depth <= 10:
            if index == 0 and depth > 0.05:
                h2_thickness += df.loc[index + 1][depth_column_name] - depth
            elif index == 0 and depth <= 0.05:
                h2_thickness += depth
            else:
                h2_thickness += depth - df.loc[index - 1][depth_column_name]
        if depth > 10:
            break


    h1_column_name = "h1_cumulative" + FS_column_name.lstrip("FS")
    h2_columnn_name = "h2_cumulative" + FS_column_name.lstrip("FS")
    df.at[0, h1_column_name] = h1_thickness
    df.at[0, h2_columnn_name] = h2_thickness

    return df

def LPI(df,depth_column_name, FS_column_name,date):
  def Integrate_LPI(z):
    return(1 - row[FS_column_name]) * (10 - 0.5 * z)

  LPI = 0
  for i, row in df.iterrows():
    depth = row[depth_column_name]
    FS = row[FS_column_name]
    if depth <=20 and FS <= 1:
      if i == 0:
        thick = df.loc[i + 1][depth_column_name] - depth
        start_depth = depth - thick
        LPI += integrate.quad(Integrate_LPI, start_depth, depth)[0]
      else:
          depth_before = df.loc[i - 1][depth_column_name]
          LPI += integrate.quad(Integrate_LPI, depth_before,depth)[0]
  df.at[0,"LPI_"+date] = LPI

  return df
def LPIish (df,depth_column_name, FS_column_name,date,h1_column_name):
  def Integrate_LPIish(z):
    c = 0
    # print(row[depth_column_name], 5/(25.56*(1-row[FS_column_name])))
    mFS = np.exp(5/(25.56*(1-row[FS_column_name])))-1
    if row[FS_column_name] <= 1 and (h1 * mFS) <= 3:
      c=(1-row[FS_column_name])
    return (25.56/z)*c
  LPIish = 0
  for i, row in df.iterrows():
    h1 = df.loc[0][h1_column_name]
    depth = row[depth_column_name]
    if h1 <= depth <= 20:
      if depth < 0.4:
        continue
      if i == 0:
        thick = df.loc[i+1][depth_column_name] - depth
        start_depth = depth - thick
        LPIish += integrate.quad(Integrate_LPIish, start_depth,depth)[0]
      else:
        LPIish += integrate.quad(Integrate_LPIish, df.loc[i-1][depth_column_name], depth)[0]
  df.at[0,"LPIish_" + date] = LPIish

  return df

def LSN(df, depth_column_name, qc1ncs_column_name, FS_column_name, date):
    def A1(qc1ncs):
        return 102 * qc1ncs ** -.82

    def A3(qc1ncs):
        return 2411 * qc1ncs ** -1.45

    def A5(qc1ncs):
        return 1701 * qc1ncs ** -1.42

    def A7(qc1ncs):
        return 1690 * qc1ncs ** -1.46

    def A9(qc1ncs):
        return 1430 * qc1ncs ** -1.48

    def A10(qc1ncs):
        return 64 * qc1ncs ** -.93

    def A11(qc1nc):
        return 11 * qc1ncs ** -.65

    def A12(qc1nc):
        return 9.7 * qc1ncs ** -.69

    def A13(qc1nc):
        return 7.6 * qc1ncs ** -.71

    def A14(qc1ncs):
        return 0

    def interpolator(lower_limit_FS_ev, lower_limit_FS, upper_limit_FS_ev, upper_limit_FS, FS):

        range = lower_limit_FS_ev - upper_limit_FS_ev

        return (upper_limit_FS - FS) * 10 * range + upper_limit_FS_ev

    def Integrate_LSN(z):
        return eps * 10 / z

    below_range_counter = 0
    LSN = 0
    total_rows_qc1ncs_qualifies = 0


    for i, row in df.iterrows():
        qc1ncs = row[qc1ncs_column_name]
        eps = np.nan
        FS = row[FS_column_name]
        depth = row[depth_column_name]
        if row[depth_column_name] <= 20 and not np.isnan(qc1ncs):
            total_rows_qc1ncs_qualifies += 1

            if 20 <= qc1ncs <= 200:
                if qc1ncs < 33:
                    eps = 10
                    below_range_counter += 1
                else:
                    eps = A1(qc1ncs)


            if .5 <= FS <= .6 and 147 <= qc1ncs <= 200:
                eps = interpolator(A1(qc1ncs), .5, A3(qc1ncs), .6, FS)

            elif .6 <= FS <= .7 and 110 <= qc1ncs <= 200:
                if qc1ncs < 147:
                    eps = interpolator(eps, .6, A5(qc1ncs), .7, FS)
                else:
                    eps = interpolator(A3(qc1ncs), .6, A5(qc1ncs), .7, FS)

            elif .7 <= FS <= .8 and 80 <= qc1ncs <= 200:
                if qc1ncs < 110:
                    eps = interpolator(eps, .7, A7(qc1ncs), .8, FS)
                else:
                    eps = interpolator(A5(qc1ncs), .7, A7(qc1ncs), .8, FS)

            elif .8 <= FS <= .9 and 60 <= qc1ncs <= 200:
                if qc1ncs < 80:
                    eps = interpolator(eps, .8, A9(qc1ncs), .9, FS)
                else:
                    eps = interpolator(A7(qc1ncs), .8, A9(qc1ncs), .9, FS)

            elif .9 <= FS <= 1 and 20 <= qc1ncs <= 200:  # check out these bounds and the ones below
                if qc1ncs < 60:
                    eps = interpolator(eps, .9, A10(qc1ncs), 1, FS)
                else:
                    eps = interpolator(A9(qc1ncs), .9, A10(qc1ncs), 1, FS)

            elif 1 <= FS <= 1.1 and 20 <= qc1ncs <= 200:
                eps = interpolator(A10(qc1ncs), 1, A11(qc1ncs), 1.1, FS)

            elif 1.1 <= FS <= 1.2 and 20 <= qc1ncs <= 200:
                eps = interpolator(A11(qc1ncs), 1.1, A12(qc1ncs), 1.2, FS)

            elif 1.2 <= FS <= 1.3 and 20 <= qc1ncs <= 200:
                eps = interpolator(A12(qc1ncs), 1.2, A13(qc1ncs), 1.3, FS)

            elif FS >= 1.3 and 20 <= qc1ncs <= 200:
                if FS > 2:
                    FS = 2
                eps = interpolator(A13(qc1ncs), 1.3, A14(qc1ncs), 2, FS)

        if i == 0:
            LSN = 0
        elif not np.isnan(eps):
            LSN += integrate.quad(Integrate_LSN, df.loc[i - 1][depth_column_name], depth)[0]

    df.at[0, "LSN_" + date] = LSN


    # print("Number of qc1ncs values below range:" ,below_range_counter,"/",total_rows_qc1ncs_qualifies)

    return df

def preforo_check(df, GWT_column_name, preforo_column_name):
    GWT_val = df.loc[0][GWT_column_name]
    preforo_val = df.loc[0][preforo_column_name]
    if GWT_val >= preforo_val:
        preforo_check = "GWT is deeper than preforo"
    elif pd.isna(preforo_val):
        preforo_check = "Nan preforo"
    else:
        preforo_check = "GWT is above preforo"
    return preforo_check
