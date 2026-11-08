# Graph-RAG Repository Standards

本仓库是五类生物数据的对外入口。首发小包与历史大包的验收状态分开记录，已有大型 Release 可继续存放在 protian-entity，由目录给出明确链接。

## 收录
- 五类结果及其实际状态：protein / rna / molecule / interaction / disease
- 对应 manifest / checksum / validation / QA 报告
- 可复现的 pipelines、contracts、脚本
- 共享 schema 和 release 工具

## 不收录
- 原始下载数据和中间产物（这些留在 protian-entity 开发仓库）
- 过程协作文档、复盘草稿、中间缓存

## 分发规则
- 稳定且中等体量：可直接入 Git
- 超大资产：优先 GitHub Release（本仓只保留元数据）
- 首发试用包：独立 `starter-v1.tar.gz` 和单行整包摘要 `SHA256SUMS.txt`；生成物进入忽略的 `build/`，不提交重复数据表
- 数据沿用上游许可与署名，代码 MIT 不覆盖数据；字段来源不明时不进入首发包

## 发布前检查
1. 明确发布对象、输入版本、字段来源、许可和筛选范围
2. `tests/release` 与对应对象的完整性检查通过；空选择、零检查和缺摘要不能作为成功
3. 首发包在干净目录校验、解包并运行三个例子；候选上传后再次匿名下载验收，再标正式
4. 修改已有域时复核其 current/index、路径和兼容调用方；不把首发小包通过写成五域历史链路全绿，也不绕过仓库实际保护规则
5. 公布的坏版本保持候选/不可用状态，修复另发新版本，不替换已公开字节

现有非蛋白 current/path 适配和旧矩阵问题保留在已知缺口，不删除 CI job 或放宽断言掩盖。新版本不隐式迁移旧发布，目录链接也不替代资产校验。
