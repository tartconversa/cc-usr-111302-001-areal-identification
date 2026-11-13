# 区划变更下的区域统计识别范围

本项目保存虚构数据共享及产品运营的合成表格与Python读取工具。`src/area_reader.py` 只解析文件；`src/areal_analysis.py` 在此之上做**可行性检查、整数端点范围与完整交叉格见证**。所有地区、机构、产品与人数均为合成，容量上界只是合成的外部公开约束。分析器不读取个人级记录或定位资料。

## 字段与单位

regions固定旧区/新区身份及总量；cells每行是允许交叉格，capacity为空表示仅受总量约束。面积权重是基准假设，不等于实际整数人数。

material-index.json记录全部输入内容哈希；fixtures保存原始资料与配置，计算输出可放.workspace-output，原始资料不应被覆盖。

分析器验证以下输入契约，冲突时直接报错而不是静默修正：

- `regions.version` / `regions.unit` 与交叉表、crosswalk 中显式声明的版本、单位必须一致；
- 旧区总量、容量、公开小计必须是非负整数；容量留空表示未设上界（与 `0` 不同）；
- cells 只允许 `old,new,capacity,area_weight`（可另有声明一致的版本/单位列）；出现 `person_id`、经纬度等疑似个人级/定位列或未知列一律拒绝；
- 每个旧区与新区都必须至少有一个允许交叉格；交叉格不得重复。

## 约束语义

- 旧区行和等于已知总数；交叉格为非负整数；**不把任何格固定为 0**。
- 容量 `<=` 上界与小计 `=` 都是硬约束。
- 公开小计是**另加的一行**，不会重复计入旧区总量，也不会被均分回各格；仅当它与总量/容量共同缩小时才改变可行域。
- 面积权重只用于独立的基准摊分（最大余数法），绝不会用来收紧、投影或解释可行范围。基准若违反容量或小计，报告 `noncompliant` 并原样展示，绝不投影后再称为面积摊分。
- 范围是**已识别的可行整数范围**，不是抽样置信区间。

## 范围、端点与见证

对每个目标区域和区域组表达式（如 `X-Y`），求解器直接优化联合表达式，返回整数最小/最大值。每个端点都附一张覆盖全部交叉格的完整见证表，并经独立校验。上下端点不会用各自边际范围拼接：`W-Y` 这类差值在边际上可能是 `[-9,0]`，联合可行范围却只有 `[-9,-1]`；出现这种差异时报告会列出 `naive_marginal_bounds` 与 `joint_bounds`。

整数求解是精确分支定界（有理数单纯形 + 精确剪枝）。若节点预算内不能认证整数端点，则报告精确的连续松弛界并保留 `unverified`，不会把松弛值冒充认证的整数端点；`solve_integer=false` 时整体标记为 `linear_relaxation_range`。无可行分配时报告删除过滤得到的冲突核（含约束来源），且**始终保留全部旧区总量**，不通过放松总量来伪造冲突。

## 运行分析器

Python 3.11以上，标准库即可，不需安装第三方依赖。

    python -m src.area_reader fixtures/regions.json
    python -m src.areal_analysis
    python -m src.areal_analysis --output .workspace-output/report.json
    python -m src.areal_analysis --no-area-baseline
    python -m compileall -q src tests
    python -m unittest discover -s tests -v

fixture 示例（A总6、B总4，A→X/Y，B→Y，A→Y容量2，公开X+Y=10）的期望结果为：

| 目标 | 最小值 | 最大值 | 说明 |
| --- | --- | --- | --- |
| X | 4 | 6 | 公开小计不把 X 缩成 5 |
| Y | 4 | 6 | 与 X 互补，但范围仍由联合约束给出 |
| X−Y | −2 | 2 | 两端各有一张可行交叉格表 |

该 fixture 的声明面积权重（A 的 1/2、1/2）给出 A→Y=3，超过容量 2，因此面积基准被标为不合规且不做投影。

## 固定容器

    docker build -t areal-identification:local .
    docker run --rm areal-identification:local

容器镜像以完整digest固定，默认运行输入契约与可行性/端点测试；实现新增分析功能所需依赖由开发者明确登记（当前分析器仅用标准库）。
