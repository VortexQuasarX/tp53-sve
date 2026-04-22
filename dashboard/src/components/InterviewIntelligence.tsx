import React from 'react'

export const InterviewIntelligence: React.FC = () => {
  const data = {
    "candidate_id": "CAND-001",
    "job_role": "Senior Python Developer",
    "scores": {
      "technical_knowledge": 8.5,
      "communication_skills": 7.0,
      "problem_solving": 8.0,
      "role_fit": 8.5,
      "confidence_sentiment": 8.0
    },
    "overall_score": 8.05,
    "recommendation": "Hire",
    "strengths": [
      "Technical expertise in Python and AWS",
      "Experience with microservices",
      "Scalability solutions using Kubernetes and Redis"
    ],
    "weaknesses": [
      "Limited detail in communication",
      "Lack of specific examples or elaboration"
    ],
    "summary": "The candidate demonstrates strong technical knowledge in Python and AWS, with relevant experience in scaling microservices using Kubernetes and Redis."
  }

  return (
    <div className="interview-module">
      <div className="module-header">
        <h1>Interview Intelligence</h1>
        <p>AI-Generated evaluation of candidate performance and technical fit.</p>
      </div>

      <div className="g2">
        <div className="card rating-summary">
          <div className="role-chip">{data.job_role}</div>
          <h2>{data.candidate_id} Evaluation</h2>
          
          <div className="score-main">
            <span className="sc-val">{data.overall_score}</span>
            <span className="sc-max">/10</span>
          </div>
          
          <div className={`rec-badge ${data.recommendation.toLowerCase()}`}>
            RECOMMENDATION: {data.recommendation.toUpperCase()}
          </div>

          <div className="score-bars">
            {Object.entries(data.scores).map(([key, val]) => (
              <div key={key} className="score-bar-row">
                <div className="sb-label">
                  <span>{key.replace('_', ' ')}</span>
                  <span>{val}/10</span>
                </div>
                <div className="sb-track">
                  <div className="sb-fill" style={{ width: `${(val as number) * 10}%` }}></div>
                </div>
              </div>
            ))}
          </div>
        </div>

        <div className="card insights-card">
          <h3><span className="ico">💡</span> Deep Insights</h3>
          <p className="summary-text">{data.summary}</p>
          
          <div className="analysis-group">
            <h4 className="green">Strengths</h4>
            <ul>
              {data.strengths.map((s, i) => <li key={i}>{s}</li>)}
            </ul>
          </div>

          <div className="analysis-group">
            <h4 className="red">Areas for Development</h4>
            <ul>
              {data.weaknesses.map((w, i) => <li key={i}>{w}</li>)}
            </ul>
          </div>
        </div>
      </div>

      <style>{`
        .role-chip { 
          background: var(--accent-light); 
          color: var(--accent); 
          padding: 4px 12px; 
          border-radius: 20px; 
          font-size: 0.75rem; 
          font-weight: 700;
          display: inline-block;
          margin-bottom: 12px;
        }
        .score-main { text-align: center; margin: 20px 0; }
        .sc-val { font-size: 3.5rem; font-weight: 900; color: var(--text-main); }
        .sc-max { font-size: 1.2rem; color: var(--text-grey); font-weight: 600; }
        
        .rec-badge { 
          text-align: center; 
          padding: 10px; 
          border-radius: 8px; 
          font-weight: 800; 
          font-size: 0.85rem; 
          margin-bottom: 24px;
        }
        .rec-badge.hire { background: #10b981; color: white; }
        
        .score-bars { display: flex; flex-direction: column; gap: 14px; }
        .sb-label { display: flex; justify-content: space-between; font-size: 0.75rem; font-weight: 600; color: var(--text-grey); text-transform: capitalize; margin-bottom: 4px; }
        .sb-track { height: 6px; background: #f3f4f6; border-radius: 3px; overflow: hidden; }
        .sb-fill { height: 100%; background: var(--accent); border-radius: 3px; }

        .summary-text { font-size: 0.9rem; color: var(--text-main); line-height: 1.6; margin-bottom: 24px; font-style: italic; }
        .analysis-group h4 { font-size: 0.8rem; text-transform: uppercase; letter-spacing: 1px; margin-bottom: 12px; padding-bottom: 4px; border-bottom: 1px solid var(--border-light); }
        .analysis-group.green h4 { color: #10b981; }
        .analysis-group.red h4 { color: #ef4444; }
        .analysis-group ul { list-style: none; padding: 0; margin-bottom: 24px; }
        .analysis-group li { font-size: 0.85rem; color: var(--text-grey); margin-bottom: 8px; display: flex; gap: 10px; }
        .analysis-group li::before { content: '→'; color: var(--text3); }
      `}</style>
    </div>
  )
}
