# 시뮬레이션 결과 고지 (Simulation notice)

이 폴더의 결과는 **합성 데이터로 생성한 교육용 결과**입니다: Python 표준 라이브러리만으로 자체 구현한 모사 코드(실제 해당 도구를 사용하지 않음).
수치(정확도, p-value, 상관계수 등)는 실제 코호트나 실제 라이브러리의 성능이 아니며 임상·연구 결론에 사용할 수 없습니다.
(Results here are produced from synthetic data for teaching; they do not reflect real-cohort results or real-tool performance.)
재현 환경: 저장소 루트의 `requirements-pipelines.txt` 참조.

**요약 지표 주의 (SCI-03-042)**: `pipeline_summary_metrics.json`의 top_sender는 out-degree 기준(myCAF 11.5168), top_receiver는 in-degree 기준(Endothelial 7.3383)입니다(2026-10-20 재생성; 이전에는 중심성 표 첫 행(out+in 합 1위, TAM_M2)을 top sender로 기록했음). 합성 데이터의 교육용 결과입니다.
