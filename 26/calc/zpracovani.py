from protokol import *

def make_df():
    global df

    df = pd.DataFrame(
        arr([
            [
            uf(38.5,  1),
            uf(53.3,  1),
            uf(77.9,  1.1),
            uf(97.9,  1.1),
            uf(112.5, 1.1),
            uf(126.2, 1.2),
            ],
            [
            uf(54.1,  1),
            uf(109.9, 1.1),
            uf(216,   1.4),
            uf(332,   1.9),
            uf(436,   2.3),
            uf(549,   2.9),
            ]
        ]).T, 
        index=[1, 2, 4, 6, 8, 10], 
        columns=["ch3cooh", "hcl"]
    )
    
    df.ch3cooh *= 1e-4
    df.hcl *= 1e-4

def add_concentrations():
    global df

    df["c_hcl"] = arr([0.0005, 0.001, 0.002, 0.003, 0.004, 0.005]) * 1e3
    df["c_ch3cooh"] = arr([0.0001, 0.0002, 0.0004, 0.0006, 0.0008, 0.001]) * 1e3

def add_molar_conductivity():
    global df

    df["sm_hcl"] = df.hcl / df.c_hcl
    df["sm_ch3cooh"] = df.ch3cooh / df.c_ch3cooh

def main():
    make_df()
    add_concentrations()
    add_molar_conductivity()
    print(df)

    chyba_mereni = arr([109.9, 111.8, 112, 110.8]).std(ddof=1)
    print(chyba_mereni)

    

if __name__ == '__main__':
    main()