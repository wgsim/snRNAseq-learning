# 시뮬레이션 결과 고지 (Simulation notice)

이 폴더의 결과는 **합성 데이터로 생성한 교육용 결과**입니다: numpy/scipy 기반 합성 데이터 시뮬레이션(실제 도구 호출 없음 또는 일부 호출).
수치(정확도, p-value, 상관계수 등)는 실제 코호트나 실제 라이브러리의 성능이 아니며 임상·연구 결론에 사용할 수 없습니다.
(Results here are produced from synthetic data for teaching; they do not reflect real-cohort results or real-tool performance.)
재현 환경: 저장소 루트의 `requirements-pipelines.txt` 참조.

**Moran's I p-값 주의 (SCI-02-067)**: `spatial_variable_genes_moran.csv`의 p-값/FDR은 순열 검정이 아니라 정규성 가정 하의 정확 분산(Cliff & Ord)을 쓴 해석적 Z 검정입니다(2026-10-20 재생성: 이전의 1/(k·n) 근사는 분산을 과소 추정해 Z를 약 1.38배 부풀렸음). 합성 데이터이며 공간 구조가 강해 p=0.0이 많으므로 유의성 수치는 해석하지 말고 I 값의 상대 순위만 참고하세요. 디콘볼루션 "벤치마크"는 참 시그니처를 solver에 제공한 합성 데이터 sanity check이며(Spearman은 희소 유형에서 0.41~0.61로 낮음), 절대 세포 수는 총 UMI에서 재스케일한 값입니다.
