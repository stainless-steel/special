mod floatext;
pub use self::floatext::FloatExt;

#[cfg(feature = "system_math")]
mod system;

#[cfg(feature = "system_math")]
pub use self::system::*;

#[cfg(not(feature = "system_math"))]
mod rust;

#[cfg(not(feature = "system_math"))]
pub use self::rust::*;
