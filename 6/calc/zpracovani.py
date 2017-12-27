from protokol import *


def init_u1():
    u1a = pd.read_csv("../raw/u1a.txt", index_col=0, sep=r"\s*,\s*")
    
    
    u1d = pd.read_csv("../raw/u1d.txt", index_col=0, sep=r"\s*,\s*")
    return u1a, u1d


def init_u4():
    p = pd.read_csv("../raw/u4p.txt", sep=r"\s*,\s*")
    s = pd.read_csv("../raw/u4s.txt", sep=r"\s*,\s*")
    return p, s


def init_u5():
    df = pd.read_csv("../raw/u5.txt", sep=r"\s*,\s*")
    return df


def main():
    u1a, u1d = init_u1()
    u4p, u4s = init_u4()
    u5 = init_u5()

    print(u1a)
    print(u1d)
    print(u4p)
    print(u4s)
    print(u5)


if __name__ == '__main__':
    main()
