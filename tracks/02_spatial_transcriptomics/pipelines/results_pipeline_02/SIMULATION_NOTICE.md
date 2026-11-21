# 시뮬레이션 결과 고지 (Simulation notice)

이 폴더의 결과는 **합성 데이터로 생성한 교육용 결과**입니다: numpy/scipy 기반 합성 데이터 시뮬레이션(실제 도구 호출 없음 또는 일부 호출).
수치(정확도, p-value, 상관계수 등)는 실제 코호트나 실제 라이브러리의 성능이 아니며 임상·연구 결론에 사용할 수 없습니다.
(Results here are produced from synthetic data for teaching; they do not reflect real-cohort results or real-tool performance.)
재현 환경: 저장소 루트의 `requirements-pipelines.txt` 참조.

**니치 해석 주의 (SCI-02-068)**: `niche_composition_summary.csv`의 니치 라벨은 데이터 기반(`Niche_<k> (<우세 세포유형>-dominant)`)이며 생물학적 이름(Cancer Core 등)은 가설로만 취급해 결과에 기록하지 않습니다(2026-10-20 재생성; 이전에는 우세 세포유형으로 생물학적 이름을 붙이는 순환 라벨링이었음). 세포 구성·수치는 동일합니다. 합성 코호트는 세포 유형 5종뿐(B 세포·DC·파골세포 없음, 골 전선에 종양 없음)이라 07·08장에 서술한 4대 니치/파골세포 활성화를 재현하지 않습니다. RANKL 축 상위 흐름도 Tumor→Tumor(79.4%)이며 "파골세포 활성화" 결과가 아닙니다. 파이프라인은 코드 흐름 시연용 합성 데모입니다.
