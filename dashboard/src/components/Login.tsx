import React from 'react'

interface LoginProps {
  onLogin: (user: any) => void
}

export const Login: React.FC<LoginProps> = ({ onLogin }) => {
  const handleGoogleLogin = () => {
    // Simulated Google Login
    const mockUser = {
      name: 'Structural Biologist',
      email: 'lab@tp53-sve.res',
      avatar: 'https://ui-avatars.com/api/?name=Strucure+Bio&background=0284c7&color=fff'
    }
    onLogin(mockUser)
  }

  return (
    <div className="login-screen">
      <div className="login-card">
        <div className="logo-icon large">p53</div>
        <h1>TP53-SVE Portal</h1>
        <p>Structural Variation Engine for Pathogenicity Prediction</p>
        
        <button className="google-btn" onClick={handleGoogleLogin}>
          <img src="https://www.gstatic.com/firebasejs/ui/2.0.0/images/auth/google.svg" alt="Google" />
          <span>Sign in with Google</span>
        </button>

        <div className="login-footer">
          <small>Secure, Encrypted Access for Enterprise</small>
        </div>
      </div>

      <style>{`
        .login-screen {
          height: 100vh;
          display: grid;
          place-items: center;
          background: var(--bg);
          position: relative;
          overflow: hidden;
        }
        .login-screen::after {
          content: '';
          position: absolute;
          width: 50vw;
          height: 50vw;
          background: radial-gradient(circle, rgba(79, 70, 229, 0.1) 0%, transparent 70%);
          top: -10vw;
          right: -10vw;
          z-index: 0;
        }
        .login-card {
          background: var(--bg2);
          border: 1px solid var(--border);
          border-radius: 24px;
          padding: 3rem;
          width: 100%;
          max-width: 440px;
          text-align: center;
          backdrop-filter: blur(20px);
          box-shadow: 0 20px 50px rgba(0,0,0,0.3);
          z-index: 1;
          animation: slideUp 0.6s ease-out;
        }
        .logo-icon.large {
          width: 80px;
          height: 80px;
          margin: 0 auto 1.5rem;
          font-size: 1.5rem;
          border-radius: 20px;
        }
        .login-card h1 {
          font-size: 2.2rem;
          margin-bottom: 0.5rem;
          background: linear-gradient(135deg, #fff, var(--accent2));
          -webkit-background-clip: text;
          -webkit-text-fill-color: transparent;
        }
        .login-card p {
          color: var(--text2);
          margin-bottom: 2.5rem;
          font-size: 0.95rem;
        }
        .google-btn {
          width: 100%;
          padding: 0.8rem;
          background: white;
          color: #333;
          border: none;
          border-radius: 12px;
          display: flex;
          align-items: center;
          justify-content: center;
          gap: 12px;
          cursor: pointer;
          font-weight: 600;
          font-size: 1rem;
          transition: all 0.2s;
        }
        .google-btn:hover {
          transform: translateY(-2px);
          box-shadow: 0 5px 15px rgba(255,255,255,0.1);
        }
        .login-footer {
          margin-top: 2rem;
          color: var(--text3);
        }
        @keyframes slideUp {
          from { opacity: 0; transform: translateY(20px); }
          to { opacity: 1; transform: translateY(0); }
        }
      `}</style>
    </div>
  )
}
