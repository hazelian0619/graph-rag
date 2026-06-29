# Graph-RAG Repository Standards

本仓库是四个生物实体域的唯一对外交付仓库。

## 收录
- 四域对外可下载结果表（按功能分区）：protein / rna / molecule / interaction
- 对应 manifest / checksum / validation / QA 报告
- 可复现的 pipelines、contracts、脚本
- 共享 schema 和 release 工具

## 不收录
- 原始下载数据和中间产物（这些留在 protian-entity 开发仓库）
- 过程协作文档、复盘草稿、中间缓存

## 分发规则
- 稳定且中等体量：可直接入 Git
- 超大资产：优先 GitHub Release（本仓只保留元数据）

## 发布前检查
1. release index 校验通过
2. consistency 校验通过
3. tests/release 全绿
4. product pointer (current.json) 的 download_url 指向本仓库
