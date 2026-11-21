# 시뮬레이션 결과 고지 (Simulation notice)

이 폴더의 결과는 **합성 데이터로 생성한 교육용 결과**입니다: Python 표준 라이브러리만으로 자체 구현한 모사 코드(실제 해당 도구를 사용하지 않음).
수치(정확도, p-value, 상관계수 등)는 실제 코호트나 실제 라이브러리의 성능이 아니며 임상·연구 결론에 사용할 수 없습니다.
(Results here are produced from synthetic data for teaching; they do not reflect real-cohort results or real-tool performance.)
재현 환경: 저장소 루트의 `requirements-pipelines.txt` 참조.

**벤치마크 수치 주의 (ENG-034)**: 이 폴더의 시간·처리량·가속배율(예: "…x speedup") 값은 단일 실행의 벽시계 측정이라 실행마다 크게 달라지며(재실행 시 673.6x → 1687.5x 관찰), 타이머 해상도 한계(예: 0.00 ms)를 분모로 한 비율은 의미가 없습니다. 배율을 인용하지 말고 "시스템 의존 예시"로만 해석하세요. (반복 측정·중앙값으로 바꾸면 커밋된 결과가 달라지므로 코드는 변경하지 않았습니다.)
