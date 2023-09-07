_sys = """假设你是一个生物医药方向的数据清洗专家，需要对很多微生物的菌名和株名进行归一化处理，请遵循进一步指导尽力提供协助。
"""

_human = """请将'{input}'里每组数据拆成菌名和株名，如果菌名存在错误，请纠正；如果你不知道如何拆解，说“不知道”。请参考下面例子输出最终结果，不要输出任何额外信息。
输入：'Sporosarcina psychrophila DSM 3; Helicobacter pylori NCTC 11637 = CCUG  17874 = ATCC 43504 = JCM 12093; Deinococus radiodurans R1 = ATCC 13939 = DSM 20539; Deinococus radioduns R1 = ATCC 13939 = DSM 20539'
输出：
Sporosarcina psychrophila, DSM 3
Helicobacter pylori, NCTC 11637 = CCUG  17874 = ATCC 43504 = JCM 12093
Deinococcus radiodurans, R1 = ATCC 13939 = DSM 20539
Deinococcus radiodurans, R1 = ATCC 13939 = DSM 20539
"""

