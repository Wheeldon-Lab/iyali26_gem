# model.xml 查询速查笔记

用 COBRApy 在命令行查 iYali26 模型里的反应 / 代谢物 / 基因。

## 0. 准备(每次先做)

```bash
cd "/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem"
source .venv/bin/activate
```

> 加载时 cobra 会打一堆 `Could not identify external compartment...` 警告,**无视它们**。
> 想屏蔽,命令末尾加 `2>/dev/null`,例如 `python3 -c "..." 2>/dev/null`。

---

## 1. 查一个反应的完整内容(最常用)

```bash
python3 -c "
import cobra
m = cobra.io.read_sbml_model('model.xml')
r = m.reactions.get_by_id('xPOOL_AC_EM')
print('name:', r.name)
print('equation:', r.reaction)
print('bounds:', r.bounds)
for mt, c in r.metabolites.items():
    print(f'  {c:+g}  {mt.id}  [{mt.name}]  formula={mt.formula}')
"
```

把 `'xPOOL_AC_EM'` 换成任何反应 ID。

---

## 2. 不知道反应 ID —— 按名字搜

```bash
python3 -c "
import cobra
m = cobra.io.read_sbml_model('model.xml')
for r in m.reactions:
    if 'pool' in (r.name or '').lower() or 'pool' in r.id.lower():
        print(r.id, '|', r.name)
"
```

把 `'pool'` 换成要搜的词(`acyl` / `biomass` / `transport` ...)。

---

## 3. 某个代谢物参与哪些反应

```bash
python3 -c "
import cobra
m = cobra.io.read_sbml_model('model.xml')
x = m.metabolites.get_by_id('m567[C_em]')
print(x.name, x.formula, x.charge)
for r in x.reactions:
    print(' ', r.id, '|', r.name, '|', r.reaction)
"
```

---

## 4. 查一个代谢物本身(化学式 / 电荷 / 注释)

```bash
python3 -c "
import cobra
m = cobra.io.read_sbml_model('model.xml')
x = m.metabolites.get_by_id('m567[C_em]')
print('name:', x.name); print('formula:', x.formula); print('charge:', x.charge)
print('annotation:', dict(x.annotation))
"
```

---

## 5. 按名字搜代谢物

```bash
python3 -c "
import cobra
m = cobra.io.read_sbml_model('model.xml')
for x in m.metabolites:
    if 'acyl-coa' in (x.name or '').lower():
        print(x.id, '|', x.name, '|', x.formula or '(空)')
"
```

---

## 6. 查一个基因(注释 + 关联反应)

```bash
python3 -c "
import cobra
m = cobra.io.read_sbml_model('model.xml')
g = m.genes.get_by_id('YALI1B29895g')
print('annotation:', dict(g.annotation))
for r in g.reactions:
    print(' ', r.id, '|', r.name, '| GPR:', r.gene_reaction_rule)
"
```

---

## 7. 查 pool / 池反应(脂质泛化的"配方")

泛化代谢物(如 `acyl-CoA_` = m567)由 pool 反应定义 —— 把若干确定链长的真实分子按生理比例混成一个池。查 acyl-CoA 池:

```bash
python3 -c "
import cobra
m = cobra.io.read_sbml_model('model.xml')
r = m.reactions.get_by_id('xPOOL_AC_EM')
print('name:', r.name)
print('equation:', r.reaction)
for mt, c in sorted(r.metabolites.items(), key=lambda kv: kv[1]):
    print(f'  {c:+.4f}  {mt.id:12} [{mt.name}]  formula={mt.formula or \"empty\"}')
"
```

列出所有 pool 反应:

```bash
python3 -c "
import cobra
m = cobra.io.read_sbml_model('model.xml')
for r in m.reactions:
    if 'pool' in (r.name or '').lower() or r.id.startswith('xPOOL'):
        print(r.id, '|', r.name)
"
```

---

## 8. 追踪一条通路(把一串反应渲染成可读链条)

### 8a. 按通路关键词搜相关反应
```bash
python3 -c "
import cobra
m = cobra.io.read_sbml_model('model.xml')
keywords = ['acyltransferase','phosphatidate','diacylglycerol','cdp-diacyl','phosphatidyl','glycerol-3-phosphate o']
for r in m.reactions:
    nm = (r.name or '').lower()
    if any(k in nm for k in keywords) and not r.boundary:
        print(r.id, '|', r.name)
"
```
(换 keywords 列表搜任意通路)

### 8b. 把一串反应按代谢物名渲染成通路(最有用)
```bash
python3 -c "
import cobra
m = cobra.io.read_sbml_model('model.xml')
chain = ['R352','R1843','R1772','R196']   # 磷脂合成主干 (Kennedy)
for rid in chain:
    r = m.reactions.get_by_id(rid)
    parts = []
    for mt, c in sorted(r.metabolites.items(), key=lambda kv: kv[1]):
        nm = (mt.name or mt.id).split('_')[0]
        coef = '' if abs(c)==1 else f'{abs(c):g} '
        parts.append((c, coef+nm))
    lhs = ' + '.join(p[1] for p in parts if p[0]<0)
    rhs = ' + '.join(p[1] for p in parts if p[0]>0)
    arrow = ' <=> ' if r.reversibility else ' --> '
    print(f'{r.id} [{r.name[:42]}]')
    print(f'   {lhs}{arrow}{rhs}')
    print()
"
```
把 `chain = [...]` 换成任意一组反应 ID 即可。

### 8c. 追一个代谢物的上下游(它从哪来、到哪去)
```bash
python3 -c "
import cobra
m = cobra.io.read_sbml_model('model.xml')
target = 'phosphatidate'
for x in m.metabolites:
    if target in (x.name or '').lower():
        print('===', x.id, x.name, '===')
        for r in x.reactions:
            role = '产物(生成)' if r.metabolites[x]>0 else '底物(消耗)'
            print(f'  [{role}] {r.id} | {r.reaction}')
        break
"
```

> **磷脂合成主干(Kennedy 通路)** 参考链:
> `R352` GPAT(装第1个脂肪酸)→ `R1843` AGPAT(装第2个)→ 磷脂酸 PA →
> 分两路:`R1772` DGAT(装第3个→三酰甘油/储油) 或 `R196` CDP-DAG synthase(→膜磷脂 PI/PS/PE/PC)。
> 每步都用泛化 `acyl-CoA` —— 这就是 lump 省反应的地方。

---

## 9. 查 / 对比 Chalmers 参考模型 (iYali v4.1.2)

参考模型在 `data/chalmers_iyali/iYali.xml`。
**注意:Chalmers 的 ID 是裸 ID,无区室后缀** —— 用 `m567` 而不是 `m567[C_em]`。

```bash
# 查 Chalmers 的某个代谢物
python3 -c "
import cobra
m = cobra.io.read_sbml_model('data/chalmers_iyali/iYali.xml')
x = m.metabolites.get_by_id('m567')
print(x.id, '|', x.name, '| formula=', x.formula or 'empty', '| charge=', x.charge)
for r in x.reactions:
    print(' ', r.id, '|', r.name, '|', r.reaction)
"
```

```bash
# 同 ID 代谢物双模型对比 formula+charge
python3 -c "
import cobra
ours = cobra.io.read_sbml_model('model.xml')
chal = cobra.io.read_sbml_model('data/chalmers_iyali/iYali.xml')
ours_x = ours.metabolites.get_by_id('m567[C_em]')
chal_x = chal.metabolites.get_by_id('m567')
print('ours :', ours_x.formula or 'empty', '| charge', ours_x.charge)
print('chal :', chal_x.formula or 'empty', '| charge', chal_x.charge)
"
```

> Chalmers 区室代号:`c`=cytoplasm, `m`=mitochondrion, `p`=peroxisome,
> `er`=ER, `e`=extracellular, `lp`=lipid particle 等(完整:把上面 model 换成 chalmers 路径跑 `.compartments`)。

---

## 备注

- **ID 格式**:代谢物 ID 带区室后缀,如 `m567[C_em]`;反应 ID 通常不带,如 `xPOOL_AC_EM`。
- **`get_by_id` 报 KeyError** = ID 打错了,先用搜索(命令 2 / 5)确认。
- **区室代号**:`C_cy`=cytoplasm, `C_mi`=mitochondria, `C_pe`=peroxisome,
  `C_er`=ER, `C_em`=endosomal membrane, `C_ex`=extracellular, 等
  (完整列表:`python3 -c "import cobra; print(cobra.io.read_sbml_model('model.xml').compartments)"`)。

### 例:xPOOL_AC_EM(acyl-CoA 池)的真实内容
该反应把 6 种确定链长的脂肪酰-CoA 按生理比例混成泛化的 `acyl-CoA_` 池:

```
0.515 oleoyl-CoA (C18:1) + 0.199 linoleoyl-CoA (C18:2)
+ 0.190 palmitoyl-CoA (C16:0) + 0.042 stearoyl-CoA (C18:0)
+ 0.003 lauroyl-CoA (C12:0) + 0.002 myristoyl-CoA (C14:0)
  <=> 0.951 acyl-CoA_   (formula 空 —— 加权平均池, 无整数化学式)
```
