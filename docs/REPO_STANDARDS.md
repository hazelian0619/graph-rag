# Protein Public Repository Standards

本仓库只服务 Protein 数据产品对外交付。

## 收录
- Protein 对外可下载结果表（按功能分区）
- 对应 manifest / checksum / validation / QA 报告
- 可复现的 protein pipelines、contracts、脚本

## 不收录
- RNA / Molecule / Interaction 资产
- 过程协作文档、复盘草稿、中间缓存

## 分发规则
- 稳定且中等体量：可直接入 Git
- 超大资产：优先 GitHub Release（本仓只保留元数据）

## 发布前检查
1. release index 校验通过
2. consistency 校验通过
3. tests/release 全绿
