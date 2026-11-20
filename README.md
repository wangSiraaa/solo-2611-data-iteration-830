# 稀释配液谱系与原样浓度回算仓（PREPQUANT）

用 Python 3.11 + SQLite 实现的**独立**配液账与样品定量仓。两个组件只通过
“导出产品”这一外部格式耦合，全部数值为 `fractions.Fraction` 精确有理数，
不连接任何设备，不依赖第三方求解器或统计推断。

```
prepquant/
  rational.py   # 有理数解析/格式化（拒绝 float）
  expr.py       # 线性表达式 const + Σ a_i·x_i
  lp.py         # 精确两阶段单纯形（Bland 规则）
  ledger.py     # 配液账：有序原子批、来源账本、产品版本、更正重放
  vault.py      # 定量仓：测量导入、联合线性区间约束、精确端点
  units.py      # 整单位（着色/无色标粒）分装档案：精确联合分布、盘查条件
  contracts.py  # 两阶段抽检合同：精确放行概率与预期抽数
  scenarios.py  # 合成场景生成器（含手算真值）
  demo.py       # 端到端演示脚本（写 out/）
  units_demo.py # 分装/回并/盘查/抽检的可复跑演示（写 out/）
  cli.py        # 命令行
tests/          # pytest（70 项）
```

## 运行

```bash
python3 -m prepquant.demo          # 全流程演示 + 手算核对，写 out/demo_report.json
python3 -m prepquant.units_demo    # 分装/回并/盘查/抽检演示，写 out/units_demo_report.json
python3 -m pytest -q               # 单元/CLI 测试
python3 -m prepquant.scenarios --out out/scenarios.json   # 仅生成场景
```

## 模型与表示

* 体积为正常数 mL，溶质质量为非负一次式，浓度 = 质量/体积，均精确有理。
* 每个容器保存 `mass = const + Σ a_i·x_i`（mg），**x_i 是各原样的未知浓度**。
  转移/混合按体积比例搬运一次式，蒸发只改体积不改质量，定容只加零质量溶剂，
  废弃清零。因此浓度回算是精确的，混合来源与中间移除都不会算错。
* 每个容器另有**来源账本**：`origin -> (仍贡献的体积, 仍贡献的质量一次式)`。
  转移按体积占比分摊来源；混合合并账本；蒸发按比例缩放体积（质量表达式不变）；
  废弃清空。`explain` 可单独交付来源与余额，不依赖任何测量。
* 标准液（已知 const）与原样（x_i）同为来源，可任意混合；系数为零的测量
  仍然参与约束（可直接导致矛盾）。

### 表达式 / 重放 / 求解方法的取舍及影响

| 选择 | 影响 |
| --- | --- |
| 以**质量**为线性货币而非浓度 | 每次操作都是一次线性运算，无需在中间步骤做除法；浓度只在查询时算一次 |
| 容器内系数**稠密**（全部已声明原样） | 实现简单、`rebase` 可在后续新增原样时对齐；导出时**只导出非零原样**，避免为无关样品编造上下界 |
| **重放**而非增量修改 | 更正 = 从目标步之前的状态重放后续全部步骤；失败则整批不落库。代价是更正 O(n)，但题面规模（≤80 容器/≤20 测量）可忽略 |
| 产品 id = `P<seq>-<容器>-<内容哈希>` | 来源与读数绑定到**内容**。未消费时可更正并让旧 id 消失；已消费时更正被拒绝，旧分析快照原样保留 |
| 精确两阶段单纯形 + Bland | 端点与可行性可被本题算术逐行复核；无需浮点容差。规模小，耗时 <1 s |
| 每次分析对所有测量**联合**求解 | 共享原样的不同混合不会被拆开解；多解时只报告精确区间，绝不取均值当作真值 |

## 命令

```bash
python3 -m prepquant.cli --ledger lab.db --vault lab.vault.db init
python3 -m prepquant.cli --ledger lab.db --vault lab.vault.db prep ops.json     # 有序原子批
python3 -m prepquant.cli --ledger lab.db explain M [--json]                     # 体积/质量式/来源
python3 -m prepquant.cli --ledger lab.db export M --out m.json                  # 产品导出
python3 -m prepquant.cli --ledger lab.db --vault lab.vault.db \
        measure --name r1 --product-file m.json --observed 28 --uncertainty 0    # 测量导入
python3 -m prepquant.cli --ledger lab.db --vault lab.vault.db \
        analyze a1 --measurement r1 --bound s1:0:100                            # 求解
python3 -m prepquant.cli --ledger lab.db --vault lab.vault.db show-analysis a1  # 历史结果
python3 -m prepquant.cli --ledger lab.db correct 3 fix.json --note "wrong top-up"
python3 -m prepquant.cli --ledger lab.db log | corrections | consumed | containers
```

操作 JSON 为对象列表，`op` 取 `standard / sample / transfer / mix / solvent /
evaporate / discard`。数值用 `"12"`、`"1/3"`、`"0.25"` 这类字符串（拒绝二进制
float）。示例：

```json
[
 {"op":"standard","container":"STD","volume":"10","concentration":"100"},
 {"op":"sample","container":"RAW","volume":"10","sample":"s1"},
 {"op":"transfer","source":"RAW","volume":"2","destination":"D"},
 {"op":"solvent","container":"D","volume":"8"},
 {"op":"mix","sources":[{"container":"D","volume":"5"},{"container":"STD","volume":"1"}],"destination":"M"},
 {"op":"evaporate","container":"M","to":"4"}
]
```

产品导出格式（定量仓只认这个，不读配液私有表；`product_id` 可省略）：

```json
{"product_id":"P6-M-…", "container":"M", "volume":"4", "unknowns":["s1"],
 "mass":{"const":"100","terms":{"s1":"1"}},
 "origins":[{"origin":"STD","volume":"1","mass":{"const":"100","terms":{}}},
            {"origin":"RAW","volume":"2/3","mass":{"const":"0","terms":{"s1":"1"}}}]}
```

## 规则与不变量

* 建立标准液/原样必须落在空容器上；容器身份稳定，已装液容器不能再被
  `standard`/`sample`/`transfer`/`mix` 当作新容器写入。
* 取出体积必须 ≤ 真实剩余体积；同一份液不能供两次；`mix` 中同一容器只可出现一次。
* 蒸发后溶质不变且体积必须仍为正、必须严格减小。
* 任一步失败 → 整批回滚（`apply_batch` 先全量重放验证，成功才写库）。
* 更正未消费步骤：从其之前状态重放；已消费产品 id 涉及的更正被拒绝，
  旧分析结果与表达式快照永不改写。
* 测量 `observed ± uncertainty` ⇒ 真浓度 ∈ `[max(0,o-u), o+u]`；乘以产品体积
  后变成质量上的两条精确线性行。分析对所有测量联立，并叠加各原样给定的
  非负上下界。输出：`feasible`/`contradiction`、各原样精确 `[min,max]`、
  是否唯一（所有原样 min==max 才标记唯一）。多解只给区间和端点处的其他
  原样取值，不取平均。
* 外部产品只要符合上述导出格式即可被定量仓直接导入。

## 合成场景与手算核对

`prepquant.scenarios` 生成 6 个场景，`demo.py` 逐一核对质量式与最终区间：

1. **两级稀释 + 标准加入 + 蒸发**：`B` 中 `100 + s1` mg / 4 mL，读数 28 ⟹ `s1=12`。
2. **共享原样的不同混合**：`M1`、`M3` 两条读数联合 ⟹ `s1=15`（联合而非分别解）。
3. **区间多解**：读数 `6±1` 于 `s1/5` ⟹ `s1∈[25,35]`，标记为非唯一。
4. **不相容读数**：同一产品 `6` 与 `7` ⟹ 矛盾。
5. **两个未知量唯一解**：两条独立读数 ⟹ `s1=10, s2=40`，均唯一。
6. **零系数矛盾**：只有标准液的产品读数 50，任何原样取值都无法满足 ⟹ 矛盾。

另在同一构建里验证：取液不足整批回滚、蒸发保持质量、更正后被消费分支拒绝、
未消费更正重放、重开数据库状态一致。

`out/demo_report.json` 是运行产物，`tests/` 中还包含：
`test_lp.py`（与 3 维暴力顶点枚举对比）、`test_ledger.py`、`test_vault.py`、
`test_cli.py`，以及整单位分装与抽检的 `test_units.py`、`test_contracts.py`、
`test_units_cli.py`、`test_units_demo.py`。

## 整单位分装档案（着色/无色标粒）

配液账处理连续体积与溶质质量；另设一个**离散**档案处理整单位标粒：每个体积
单位恰含一粒，每粒或**着色**或**无色**。母盒声明两类粒数，分装时从某盒等概率、
不放回地取 `k` 粒，其余留在原盒。仅支持至多 16 粒、同时至多 4 个非空盒。

### 表示与相关性

* **联合状态，而非边缘。** 从同一母盒分出的两盒负相关：第一盒着色多，第二盒
  必然少。因此状态是**所有非空盒着色数向量的精确联合分布**
  `((盒, 着色数), …) → Fraction`，而不是各盒边缘。这样上游拆分的相关性被完整保留。
* **可交换性使计数向量充分。** 盒内哪些具体粒在哪无关紧要，只看着色数。从
  `(t, c)` 取 `k` 粒得 `d` 粒着色的概率为超几何
  `C(c,d)·C(t-c,k-d)/C(t,k)`，余下 `c-d` 粒。后续在“盒当时实际状态”条件下的
  抽取，正是在当前联合分布中对 `c` 求条件，于是连续多次拆分由**逐次超几何卷积**
  精确处理——测试中以多元超几何分布的闭式与逐次拆分的联合分布逐项比对。
* **回并是确定的。** 把分出物重新合盒只是把着色数相加，不引入随机性，但合并盒
  与其余盒的协方差仍被保留。
* **盘查即条件；矛盾是一等公民。** 非破坏性盘查只把某盒着色数限定在区间
  `[lo, hi]`，把匹配的联合概率质量重新归一。若盘查后条件概率为零、或要求超过
  该盒实有粒数、或 `lo > hi`，一律以 `UnitContradiction` 明确拒绝，绝不静默丢弃
  冲突。

记录为有序、原子、可重放的 `beads / transfer / mix / inspect` 批；任何违规都会
整批回滚。`transfer`/`mix` 与配液账同形（也接受 `volume` 别名），`beads_from_ledger`
可把配液账中真实的 `transfer`/`mix` 记录（给定母盒两类粒数、体积视为整单位数）
直接驱动本档案；若某标粒盒被 `solvent`/`evaporate`/`discard` 或与非标粒容器混用，
则明确报错，避免离散档案与配液账静默偏离。

### 命令行

```bash
python3 -m prepquant.cli --units beads.db units-prep records.json
python3 -m prepquant.cli --units beads.db units-boxes --json
python3 -m prepquant.cli --units beads.db units-joint A B --out product.json   # 两盒联合产品
python3 -m prepquant.cli --units beads.db units-log
python3 -m prepquant.cli contract --product product.json --contract c.json     # 抽检
```

记录 JSON 示例（`beads / transfer / mix / inspect`）：

```json
[
 {"op":"beads","box":"P","colored":"2","colorless":"2"},
 {"op":"transfer","source":"P","units":"2","destination":"A"},
 {"op":"mix","sources":[{"box":"A"},{"box":"B"}],"destination":"M"},
 {"op":"inspect","box":"A","min":"1","max":"2"}
]
```

`units-joint A B` 导出的**联合产品**只含两盒粒数与精确联合分布，检测端只读它与
自己的合同，不接触档案私有账本。

## 两阶段抽检合同

合同字段（均为整数，满足 `a < r`）：第一阶段从第一盒不放回抽 `n` 粒，着色数 `x`
满足 `x <= a` 提前放行、`x >= r` 提前退回，其余继续；第二阶段从第二盒不放回抽
`m` 粒，两次合计 `X+Y <= b` 才最终放行。`second_box` 可指定另一盒；缺省
`"continue"` 表示**在第一盒继续抽其余量**（即用第一阶段的剩余实物，而非重抽整盒）。

对每个第一阶段着色数 `x`，报告给出该分支概率、早放行/早退回/继续的判断、以及经
该分支最终放行的概率；再给出总放行概率与预期抽取粒数（`n + m·P(继续)`），全部
为精确分数。若某条**正概率的继续分支**需要抽取超过所选第二阶段盒实有粒数，合同
被拒绝——该分支确会发生，承诺取超过实有的粒数无效。

```
$ python3 -m prepquant.cli contract --product product.json --contract c.json
stage 1: draw 1 from 'A'; x<= 0 accept, x>= 2 reject, else continue
stage 2: draw 1 from 'P'; final accept iff total <= 1
  x   P(X=x)   branch     eventual accept
  0       1/2   accept     1/2
  1       1/2   continue   1/3
stage-1 split: accept 1/2, reject 0, continue 1/2
release probability: 5/6
expected units drawn: 3/2
```

### 手算例子（`units_demo.py` 逐一核对）

* 四粒两着色，随机分 A、B 各两粒；先从 A 抽一粒，无色即放行，着色则从 B 抽一
  粒，累计至多一着色才放行（`a=0, r=2, b=1`）：联合分布
  `{(0,2):1/6,(1,1):2/3,(2,0):1/6}`，第一阶段 `x=0`（1/2）即放行，`x=1`（1/2）
  继续并再放行 2/3，故总放行 `1/2+1/3=5/6`，预期抽取 `1 + 1·(1/2)=3/2`。
* 固定分装，改 `b=0`：`x=1` 分支不再放行，总放行 `1/2`。
* 起始档改为四粒仅一着色：总放行 `1`（必放行），预期抽取 `5/4`。
* **同一盒续抽**：四粒两着色、A 两粒，先从 A 抽一粒，着色则在 A 用其仅剩的一粒
  继续抽（`second_box="continue"`）。这是对余量取样而非重抽整盒，总放行 `5/6`
  （重抽整盒会错误得到 `2/3`）。

## 精确性保证

* 输入禁止 `float`；内部只用 `Fraction`。
* 单纯形返回的每个顶点都会用原始约束逐行 `Fraction` 复核（`_verify`），
  不一致即抛错而不是给出近似值。
* `tests/test_lp.py` 在二维/三维随机系统上把可行性、每个变量的最小/最大
  端点与暴力枚举全部顶点（Cramer 法则）比对。
* 整单位档案全程精确：联合分布、超几何概率、盘查重归一、抽检放行概率与预期
  抽取粒数都是 `Fraction`；逐次拆分与多元超几何分布闭式逐项比对，抽检三个示例
  的放行概率/抽数与手算一致。
