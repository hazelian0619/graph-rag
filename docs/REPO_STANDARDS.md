# Graph-RAG Public Repository Standards

本仓库定位：**对外结果产物仓库（Result-first）**。

## 1. 收录标准（必须）
- 可复用：对外用户可直接下载、验证、消费
- 可追溯：必须能关联到 manifest / checksum / QA 报告
- 可维护：路径稳定（release/index.json + products/*/current.json）

## 2. 不收录标准（禁止）
- 过程性协作材料（复盘草稿、分工记录、会议纪要）
- 临时中间产物、缓存、实验脚本副本
- 无法复现来源或未通过验证门禁的数据

## 3. 数据分发规则
- <=100MB 且稳定核心表：可放 Git
- >100MB 或高频变更大文件：放 GitHub Release 资产
- 仓库内仅保留 Release 元数据（manifest / checksums / QA JSON）

## 4. 结构主线（四产品）
- protein
- rna
- molecule
- interaction

## 5. 对外检查清单（发布前）
1. `python3 scripts/validate_release_index.py ...` PASS
2. `python3 scripts/check_release_consistency.py ...` PASS
3. `pytest -q tests/release` PASS
